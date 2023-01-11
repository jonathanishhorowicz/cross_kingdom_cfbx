library(MOFA2)
reticulate::use_condaenv("R4")
library(stringr)
library(tidyr)
library(data.table)

source("../scripts/utils.R")

###############################################################
# load data
###############################################################

sample_groups <- readRDS("../data/sample_groups.rds")

#
# class labels for each task
y <- fread("../data/phenotypes.csv") %>%
  as_tibble() %>%
  column_to_rownames("sample_id")

# remove samples with missing data
y <- y[!rownames(y) %in% c("FAME00008", "FAME00381"),]

# phyloseqs
phylos <- readRDS("../data/final_phyloseqs.rds")
phylos <- lapply(phylos, prune_uncultured)
phylos <- lapply(phylos, fill_in_taxonomy)

###############################################################
###############################################################

###############################################################
# prepare aggloemerated OTU tables
###############################################################

TAXA_ORDER <- c("Class", "Order", "Family", "Genus", "Species")

phylos_glom <- lapply(
  setNames(TAXA_ORDER, TAXA_ORDER),
  function(rank_name) lapply(phylos, function(x) speedyseq::tax_glom(x, rank_name))
) 

# agglomerated OTU tables
X <- lapply(
  phylos_glom,
  function(x) lapply(x, function(xx) xx %>% otu_table() %>% as.data.frame() %>% t())
) %>% unlist(recursive=FALSE)
X <- lapply(X, function(x) x[ rownames(y), ])
lapply(X, dim)

# CLR transform
X_clr <- lapply(X, compositions::clr)
lapply(X_clr, dim)

# Remove rare OTUs (those with prevalence below 10%)
taxa_prev <- lapply(X, function(x) 100*apply(x, 2, function(x) sum(x>0)/length(x)))
rare_taxa <- lapply(taxa_prev, function(x) x[x < 10])
rare_taxa <- unlist(rare_taxa)
names(rare_taxa) <- names(rare_taxa) %>%
  str_remove_all("bacterial\\.|fungal\\.") %>%
  str_remove_all(paste0(TAXA_ORDER, "\\.", collapse="|"))
X_clr <- lapply(X_clr, function(x) x[ , !colnames(x) %in% names(rare_taxa)] )
lapply(X_clr, dim)

# check sample IDs match the phenotypes
sapply(X_clr, function(x) identical(rownames(x), rownames(y))) %>%
  all() %>%
  stopifnot()

###############################################################
###############################################################

###############################################################
# run MOFA on each agglomerated OTU table
###############################################################

MOFA_var_explained <- list()

for(taxa_rank in TAXA_ORDER) {
  cat(taxa_rank, "\n")
  
  names(X_clr)[ grepl(taxa_rank, names(X_clr)) ]
  
  MOFAobject <- create_mofa(lapply(X_clr[ grepl(taxa_rank, names(X_clr)) ], t)) # , groups=y$Disease)
  mofa_opts <- get_default_data_options(MOFAobject)
  mofa_opts$scale_views <- TRUE
  MOFAobject <- prepare_mofa(MOFAobject, data_options=mofa_opts)
  MOFAobject <- run_mofa(MOFAobject)
  
  MOFA_var_explained[[ taxa_rank ]] <- MOFAobject %>%
    calculate_variance_explained()
}

###############################################################
###############################################################

###############################################################
# plot the variance explained by the factors at each rank
###############################################################

# collect variance explained for the plots
mofa_var_plotdata <- lapply(
  MOFA_var_explained,
  function(x) x$r2_per_factor$group1 %>%
    as.data.frame() %>%
    rownames_to_column("Factor") %>%
    magrittr::set_colnames(str_remove(colnames(.), paste0(TAXA_ORDER, "\\.", collapse="|")))
  ) %>%
  bind_rows(.id="Rank") %>%
  pivot_longer(c(fungal, bacterial), names_to="view") %>%
  mutate(
    Rank=factor(Rank, levels=TAXA_ORDER),
    view=recode(view, fungal="Fungal\nview", bacterial="Bacterial\nview")
    )

# plot the variance explained for al the agglomeration ranks
mofa_var_plotdata %>%
  mutate(
    Factor=factor(Factor, levels=paste0("Factor", 1:12)),
    Rank=factor(Rank, levels=c("Class", "Order", "Family", "Genus", "Species"))
  ) %>%
  ggplot(aes(x=Factor, y=view, fill=value)) +
  geom_tile() +
  ggforce::facet_row(~Rank, space="free", scales="free_x") +
  scale_fill_gradient(low="gray97", high="darkblue",
                      limits=c(0, NA),
                      guide=guide_colourbar(frame.colour="black",
                                            ticks.colour="black")) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_minimal() +
  theme(
    panel.border=element_rect(fill=NA),
    axis.text.x=element_text(size=10, angle=90, hjust=1, vjust=0.5),
    axis.ticks=element_line(),
    axis.text.y=element_text(size=10),
    strip.text.x=element_text(size=12),
    axis.title.y=element_blank(),
    axis.title.x=element_text(size=14),
    legend.text=element_text(size=12),
    legend.position="right"
  ) +
  labs(fill="Var (%)") 

###############################################################
###############################################################
