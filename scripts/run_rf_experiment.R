library(dplyr)
library(data.table)
library(stringr)
library(tibble)
library(tidyr)
library(parallel)
library(foreach)
library(ranger)

options(stringsAsFactors=FALSE)

source("../scripts/utils.R")

# Read command line arguments
if(!interactive()) {
  PBS_JOBID <- Sys.getenv("PBS_JOBID")
  PBS_ARRAY_INDEX <- Sys.getenv("PBS_ARRAY_INDEX") %>% as.integer()
  N_CORES <- Sys.getenv("N_CORES") %>% as.integer()
  N_RANGER_THREADS <- Sys.getenv("N_RANGER_THREADS") %>% as.integer()
  DEBUG <- FALSE
  N_PERMS <- Sys.getenv("N_PERMS") %>% as.integer()
  VARIMP_RUN <- Sys.getenv("VARIMP_RUN")=="TRUE"
} else {
  PBS_JOBID <- "DEBUG"
  PBS_ARRAY_INDEX <- 1
  N_CORES <- 1
  DEBUG <- TRUE
  N_PERMS <- 10
  N_RANGER_THREADS <- 1
  VARIMP_RUN <- FALSE
}

# if running a PBS array job this will get the base ID so all results
# are saved in the same directory
base_jobid <- str_remove_all(PBS_JOBID, "\\.pbs|\\[\\d+\\]")

message("base_jobid: ", base_jobid)
message("PBS_ARRAY_INDEX: ", PBS_ARRAY_INDEX)
message("N_CORES: ", N_CORES)
message("N_RANGER_THREADS:", N_RANGER_THREADS)
message("VARIMP_RUN:", VARIMP_RUN)

save_path <- Sys.getenv("SAVE_PATH")
if(save_path=="") {
  save_path <- file.path("../results", base_jobid)
  cat("No save path specified - using", save_path, "\n")
}

# make save directory
dir.create(file.path(save_path))
dir.create(file.path(save_path, "settings", "job_args"), recursive=TRUE)
dir.create(file.path(save_path, "job_completion"))
dir.create(file.path(save_path, "validation_fold_preds"))
dir.create(file.path(save_path, "var_imps"))

cat(sprintf("Made save directories in %s\n", save_path))

###############################################################
# base simulation settings
################################################################

SEED <- 69012365
set.seed(SEED)

N_TREES <- ifelse(DEBUG, 100, 1000) # number of trees in random forest
N_HPARAM_SEARCH_ITER <- 10 # size of random forest hyperparamter search
N_INNER_FOLDS <- ifelse(DEBUG, 2, 5) # number of CV folds (k)

# save experiment settings
settings <- list(
  SEED=SEED,
  N_TREES=N_TREES,
  N_CORES=N_CORES,
  N_RANGER_THREADS=N_RANGER_THREADS,
  N_PERMS=N_PERMS,
  VARIMP_RUN=VARIMP_RUN,
  N_RANGER_THREADS=N_RANGER_THREADS,
  N_HPARAM_SEARCH_ITER=N_HPARAM_SEARCH_ITER,
  N_INNER_FOLDS=N_INNER_FOLDS
)

settings %>%
  as.data.frame() %>%
  fwrite(
    file.path(save_path, "settings", sprintf("%d.csv", PBS_ARRAY_INDEX))
  )

message("SEED:", SEED)
message("N_TREES:", N_TREES)
message("save_path:", save_path)
message("N_HPARAM_SEARCH_ITER:", N_HPARAM_SEARCH_ITER)
message("N_INNER_FOLDS:", N_INNER_FOLDS)

###############################################################
###############################################################

###############################################################
# this simulation settings
###############################################################

all_job_args <- expand_grid(
  phenotype=c("Disease", "Group", "Exacing"),
  kingdom=c("fungal", "bacterial", "all"),
  transform=c("identity", "clr", "log1p", "sum_norm", "relative_abundance"),
  spec_name=c("impurity_corrected--altmann", "permutation--altmann", "impurity--NA", "none--NA"),
  n_perms=c(1000),
  glom_rank=c("Class", "Order", "Family", "Genus", "Species"),
  k=c(5,10),
  min_taxa_prev=c(0,10,20)
) 

if(VARIMP_RUN) {
  # only compute variable importances when there is a detectable difference
  # between the two groups
  cat("This is a varimp run\n")
  all_job_args <- all_job_args %>%
    filter(k==5, spec_name!="none--NA", glom_rank=="Genus") %>%
    right_join(
      tibble(phenotype=c("Disease", "Disease", "Disease", "Group"),
             kingdom=c("fungal", "bacterial", "all", "fungal")),
      by=c("phenotype", "kingdom")
    )
} else {
  all_job_args <- all_job_args %>%
    filter(spec_name=="none--NA")
}

all_job_args %>% as.data.frame() %>% fwrite(
  file.path(save_path, "settings", "job_args", sprintf("%d.csv", PBS_ARRAY_INDEX))
)

# extract relevant phenotype-kingdom combinations
cat(sprintf("This is job %d of %d\n", PBS_ARRAY_INDEX, nrow(all_job_args)))

# run experiment
job_args <- all_job_args[PBS_ARRAY_INDEX,] %>% as.list()

cat("phenotype:", job_args$phenotype, "\n")
cat("kingdom:", job_args$kingdom, "\n")
cat("transform:", job_args$transform, "\n")
cat("n_perms:", job_args$n_perms, "\n")
cat("spec_name:", job_args$spec_name, "\n")
cat("resampling_type:", job_args$resampling_type, "\n")
cat("glom_rank:", job_args$glom_rank, "\n")
cat("k:", job_args$k, "\n")
cat("min_taxa_prev:", job_args$min_taxa_prev, "\n")

###############################################################
###############################################################

###############################################################
# load data
###############################################################

# samples for each of the three tasks
sample_groups <- readRDS("../data/sample_groups.rds")

# class labels for each task
y <- fread("../data/phenotypes.csv") %>%
  as_tibble() %>%
  column_to_rownames("sample_id")

#
# remove samples with missing data
y <- y[!rownames(y) %in% c("FAME00008", "FAME00381"),]

#
# load phyloseq objects (one per kingdom) and remove any uncultured taxa
phylos <- readRDS("../data/final_phyloseqs.rds")
phylos <- lapply(phylos, prune_uncultured)

# agglomerate to the required rank (after replacing any NAs in the taxonomy with
# the highest available rank)
cat(sprintf("Agglomerating to %s\n", job_args$glom_rank))
phylos <- lapply(phylos, fill_in_taxonomy)
phylos_glom <- lapply(phylos, function(x) speedyseq::tax_glom(x, job_args$glom_rank))

###############################################################
###############################################################

###############################################################
# set up random forest model fitting
###############################################################

model_specs <- list()
for(importance in c("impurity_corrected", "permutation")) {
  model_specs[[ sprintf("%s--%s", importance, "altmann") ]] <- ranger_modelspec_factory(importance, "altmann")
}
model_specs[[ "impurity--NA" ]] <- ranger_modelspec_factory("impurity", "NA")
model_specs[[ "none--NA" ]] <- ranger_modelspec_factory("none", "NA")

# some transformations (e.g. clr) needs to be applied to fungal and bacterial counts separately
if(job_args$transform %in% c("clr", "sum_norm", "relative_abundance")) {
  if(job_args$kingdom=="all") {
    split_point <- lapply(phylos_glom, function(x) x %>% otu_table() %>% nrow())[[1]]
    cat(sprintf("Splitting at column %d\n", split_point))
    transform_fn <- function(x) {split_apply(x, transforms[[ job_args$transform ]], split_point)}
  } else {
    transform_fn <- transforms[[ job_args$transform ]] 
  }
} else {
  transform_fn <- transforms[[ job_args$transform ]]
}

#
# setup design matrices from OTU tables
if(job_args$kingdom!="all") {
  included_phylo_names <- job_args$kingdom
} else {
  included_phylo_names <- c("fungal", "bacterial")
}

sample_ids <- sample_groups[[ job_args$phenotype ]]
x <- lapply(
  phylos_glom[included_phylo_names],
  function(obj) obj %>%
    otu_table(taxa_are_rows) %>%
    as.data.frame() %>%
    t()
)
x <- lapply(x, function(xx) xx[ names(sample_ids),  ])
x <- do.call(cbind, x)
cat("x has shape", nrow(x), "x", ncol(x), "\n")

# transform x as required
x_transformed <- transform_fn(x)
stopifnot(identical(colnames(x_transformed), colnames(x)))

#
# remove rare taxa (if required). rarity is computed on the 
# original X in case the transformation doesn't preserve zeroes
if(job_args$min_taxa_prev > 0) {
  cat(sprintf("Removing taxa with prevalence less than %d %%\n", job_args$min_taxa_prev))
  taxa_prev <- 100*apply(x, 2, function(x) sum(x>0)/length(x))
  rare_taxa <- taxa_prev[ taxa_prev < job_args$min_taxa_prev ]
  cat(sprintf("Removing %d taxa\n", length(rare_taxa)))
  x_transformed <- x_transformed[ ,!colnames(x_transformed) %in% names(rare_taxa) ]
}

yy <- y[ names(sample_ids), job_args$phenotype ] %>%
  as.factor() %>%
  droplevels()

# training-validation split and the permuted labels for assessing statistical significance
# Note that these are pre-calculated so we can split over different nodes on the cluster
# without worrying about repeating permutations or using the same splits
if(!VARIMP_RUN) {
  training_idxs <- sprintf("../data/kfold_cv_idxs/k=%d/%s.rds", job_args$k, job_args$phenotype) %>%
    readRDS()
  cat(length(training_idxs), "training folds\n")

  permuted_labels <- file.path("../data/permutation_idxs", sprintf("%s.rds", job_args$phenotype)) %>%
    readRDS()
  
} else {
  training_idxs <- NULL
  permuted_labels <- NULL
}

# reduce the number of permutations in the label randomisation test (for initial/debugging runs)
if(!is.na(N_PERMS) & length(permuted_labels)>0) {
  if(!DEBUG) {
    warning("Reducing number of permutations in the label randomisation test but DEBUG=FALSE",
            call.=FALSE)
  }
  permuted_labels <- permuted_labels[1:N_PERMS]
}
cat(length(permuted_labels), "sets of permuted_labels\n")

###############################################################
###############################################################

###############################################################
# main loop - fits all the required models #(observed labels
# and permuted ones for the label randomisation test)
###############################################################

iterator_list <- c("observed", paste0("permutation_", 1:length(permuted_labels)))
cat(sprintf("Running the main loop (%d iterations) with %d core(s)\n", length(iterator_list), N_CORES))

cl <- makeCluster(N_CORES, type="FORK")
registerDoParallel(cl)

dummy_var <- foreach(i=1:length(iterator_list)) %dopar%  {
  
  iteration_label <- iterator_list [[ i ]] # either "observed" or "permutationX" where X is an integer

  # permute the labels if requred
  if(iteration_label!="observed") {
    resample_idx <- iteration_label %>% 
      str_extract("\\d+") %>%
      as.integer()
    yy_sample_order <- permuted_labels[[ resample_idx ]]

    if(VARIMP_RUN) {
      x_sample_order <- permuted_labels[[ resample_idx ]]
    } else {
      x_sample_order <- 1:nrow(x)
    }

  } else{
    yy_sample_order <- 1:length(yy)
    x_sample_order <- 1:nrow(x)
  }
  
  # train a single model on a single X,y
  this_it_out <- run_experiment(
    X=x_transformed[ x_sample_order, ],
    y=yy[ yy_sample_order ],
    model_fns=model_specs[[ job_args$spec_name ]],
    training_idxs=training_idxs,
    seed=SEED,
    parallel_backend="FORK",
    n_cores=1,
    n_ranger_threads=N_RANGER_THREADS
  )
  
  # save outer fold validation predictions (if not a varimp run)
  if(!is.null(this_it_out$val_preds_labels)) {
    this_it_out$val_preds_labels %>%
      mutate(source=iteration_label) %>%
      append_sim_args(job_args) %>%
      fwrite(
        file.path(save_path, "validation_fold_preds", sprintf("%s____%s.csv", iteration_label, PBS_ARRAY_INDEX))
      )
  }
  
  # save variable importance resultss
  if(nrow(this_it_out$var_imps)>0) {
    this_it_out$var_imps %>%
      mutate(source=iteration_label) %>%
      append_sim_args(job_args) %>%
      fwrite(
        file.path(save_path, "var_imps", sprintf("%s____%s.csv", iteration_label, PBS_ARRAY_INDEX))
      )
  }
}

stopCluster(cl)
registerDoSEQ()

# log that this job of the array complete successfully
message("Script finished successfully")
file.create(file.path(save_path, "job_completion", paste0(PBS_ARRAY_INDEX, ".csv")))

###############################################################
###############################################################
