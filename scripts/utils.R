library(phyloseq)
library(caret)
library(foreach)
library(doParallel)
library(data.table)
library(dplyr)
library(tibble)

append_sim_args <- function(df, sim_args) {
  # add columns to df for each element of sim_args
  tmpdf <- matrix(rep(sim_args, nrow(df)), ncol=length(sim_args), byrow=TRUE) %>%
    as.data.frame()
  colnames(tmpdf) <- names(sim_args)
  return(df %>% bind_cols(tmpdf))
}

prune_uncultured <- function(x) {
  # Removes any OTUs whose name contains the string `uncultured`
  otu_names <- taxa_names(x)
  taxa_to_remove <- !grepl("uncultured", otu_names)
  if(length(taxa_to_remove)>0) {
    cat(sprintf("Removing %d uncultured taxa\n", length(otu_names)-sum(taxa_to_remove)))
    return(prune_taxa(taxa_to_remove, x))
  } else {
    return(x)
  }
}

fill_in_taxonomy  <- function(phylo) {
  # for a phyloseq, replace any NAs in the taxonomy table with the next available
  # (non-NA) rank
  
  new_taxonomy_mat <- tax_table(phylo)@.Data %>%
    as.data.frame() %>%
    mutate(
      Order=na_if(Order, "o__"),
      Family=na_if(Family, "f__"), Family=if_else(grepl("unidentified", Family), NA_character_, Family),
      Genus=na_if(Genus, "g__"), Genus=if_else(grepl("unidentified", Genus), NA_character_, Genus),
      Species=na_if(Species, "s__")
    ) %>% 
    mutate(
      Phylum=if_else(is.na(Phylum), Domain, Phylum),
      Class=if_else(is.na(Class), Phylum, Class),
      Order=if_else(is.na(Order), Class, Order),
      Family=if_else(is.na(Family), Order, Family),
      Genus=if_else(is.na(Genus), Family, Genus),
      Species=if_else(is.na(Species), Genus, Species)
    ) %>%
    as.matrix()
  
  tax_table(phylo) <- new_taxonomy_mat
  return(phylo)
}

# list of transformations for OTU counts
transforms <- list(
  identity=function(x) x,
  clr=function(x) as.matrix(as.data.frame(compositions::clr(x))),
  log1p=log1p,
  sum_norm=function(x) { return(x/sum(x)) },
  relative_abundance=function(x) { return(x/rowSums(x) )}
)

# returns functions to fit a random forest model and evaluate its
# variable importance
ranger_modelspec_factory <- function(
    importance=c("impurity_corrected", "permutation", "impurity", "none"),
    pvalues=c("altmann", "NA")) {
  
  importance <- match.arg(importance)
  pvalues <- match.arg(pvalues)
  
  if(importance %in% c("impurity") & pvalues!="NA") { 
    stop(sprintf("No p-values with %d importance", importance), call.=FALSE)
  }
  
  model_spec <- list(
    train=function(X, y) {
      caret::train(
        X, y,
        method="ranger",
        importance=importance,
        trControl=trainControl(method="CV",
                               number=N_INNER_FOLDS,
                               search="random",
                               classProbs=TRUE,
                               summaryFunction=twoClassSummary,
                               allowParallel=FALSE),
        metric="ROC",
        tuneLength=N_HPARAM_SEARCH_ITER,
        num.trees=N_TREES,
        num.threads=N_RANGER_THREADS,
        keep.inbag=TRUE)
    }
  )
  
  if(pvalues=="altmann") {
    model_spec$varimp <- function(fit, X, y) {
      out <- tryCatch({
        df <- cbind(X,y)
        importance_pvalues(fit$finalModel,
                           method="altmann",
                           num.permutations=job_args$n_perms,
                           formula=y~., data=df,
                           num.threads=N_RANGER_THREADS) %>%
          as.data.frame() %>%
          rownames_to_column("variable")
      }, error=function(e) {
        cat("Error in varimp")
        message(e)
        return(data.frame())
      })
    }
  } else if (pvalues=="NA") {
    if(importance=="impurity") {
      model_spec$varimp <- function(fit, X, y) {
        fit$finalModel$variable.importance %>%
          as.data.frame() %>%
          rownames_to_column("variable") %>%
          rename(importance=".")
      }
    } else if(importance=="none") {
      model_spec$varimp <- function(fit, X, y) {
        return(tibble())
      }
    }
  }
  
  return(model_spec)
}

run_experiment <- function(X, y,
                           model_fns,
                           training_idxs,
                           seed,
                           verbose=FALSE,
                           parallel_backend=c("FORK", "PSOCK"),
                           n_cores=1,
                           n_ranger_threads=1) {
  # input checks
  if(is.null(X))
    stop("X cannot be NULL")
  if(is.null(y))
    stop("y cannot be NULL")
  
  if(sum(is.na(X))!=0) stop(sprintf("X contains %d NA values", sum(is.na(X))))
  if(sum(is.na(y))!=0) stop(sprintf("y contains %d NA values", sum(is.na(y))))
  stopifnot(n_cores>0)
  
  parallel_backend <- match.arg(parallel_backend)
  
  positive_class <- levels(y)[[2]]
  if(verbose) cat(positive_class, "is the positive class\n")
  
  out <- list() # stores all output
  
  # run the cross-validation (if training indices are provided, otherwise go
  # straight to full model fitting)
  if(length(training_idxs)>0) {
    
    # Setup parallelisation and RNG seeding
    if(n_cores>1) {
      if(.Platform$OS.type=="windows" & parallel_backend=="FORK") {
        warning("FORK parallel backend has no effect on Windows. Switching to PSOCK")
      }
      if(verbose) cat(sprintf("Registering %s cluster with %d cores\n", parallel_backend, n_cores))
      cl <- parallel::makeCluster(n.cores, type=parallel_backend)
      if(parallel_backend=="PSOCK")
        clusterExport(cl, list("cv.glmnet", "factor2int", "rbindlist", "simple_roc")) # outdated
      registerDoParallel(cl)
    }
    
    if(verbose) cat("Running nested CV with", length(training_idxs), "outer folds...")
    
    cv_results <- foreach(i=1:length(training_idxs)) %dopar% {
      set.seed(seed) # make sure each loop uses the same seed - not doing any data splitting inside the loop
                     # so this only makes sure the RF model construction is reproducible
      
      # setup X and y
      X_train <- X[ training_idxs[[i]], ]
      y_train <- y[ training_idxs[[i]] ]
      X_val <- X[ -training_idxs[[i]], ]
      y_val <- y[ -training_idxs[[i]] ]
      
      stopifnot(sum(is.na(X_train))==0)
      stopifnot(sum(is.na(X_val))==0)
      stopifnot(sum(is.na(y_train))==0)
      stopifnot(sum(is.na(y_val))==0)
      
      # fit model on training set
      selected_model <- model_fns$train(X_train, y_train)
      
      # for estimating CIs of CV AUC estimate
      val_preds <- predict(selected_model,
                           newdata=X_val,
                           type=c(Regression="raw", Classification="prob")[[ selected_model$modelType ]],
                           num.threads=n_ranger_threads)
      if(selected_model$modelType=="Classification") {
        val_preds <- val_preds[,which(colnames(val_preds)==positive_class)]
        val_labels <- as.integer(y_val==positive_class)
      } else if(selected_model$modelType=="Regression") {
        val_labels <- y_val
      } else {
        stop(sprintf("Unrecognised modelType: %s", selected_model$modelType), call.=FALSE)
      }
      
      # return predictions and selected hyperparameters
      list(cv_auc_args=list(predictions=val_preds, labels=val_labels, sample_id=rownames(X_val)),
           best_tune=selected_model$bestTune)
    }
    if(verbose) cat("done\n")
    
    # close workers
    if(n_cores > 1)
      stopCluster(cl)
    registerDoSEQ()
  } else {
    if(verbose) cat("Not running any CV\n")
    cv_results <- NULL
  }
  
  # combine results from CV loops into a single dataframe
  if(length(training_idxs)>0) {
    out$best_tune <- lapply(cv_results, function(x) x$best_tune)
    names(out$best_tune) <- paste0("resample", 1:length(training_idxs))
    out$val_preds_labels <- lapply(
      cv_results,
      function(x) x$cv_auc_args %>% as_tibble()
    ) %>%
      rbindlist(idcol="outer_fold")
  } else {
    out$best_tune <- NULL
    out$val_preds_labels <- NULL
  }
  
  # then fit on the whole dataset
  stopifnot(sum(is.na(X))==0)
  
  if(verbose) cat("Fitting full dataset model...")
  full_data_model <- model_fns$train(X, y)
  if(verbose) cat("done\n")
  
  if(verbose) cat("Calculating variable importance...")
  out$var_imps[["full_dataset"]] <- model_fns$varimp(full_data_model, X, y)
  if(verbose) cat("done\n")
  
  out$best_tune[["full_dataset"]] <- full_data_model$bestTune
  
  # format outputs
  out$var_imps <- rbindlist(out$var_imps, idcol="source")
  out$best_tune <- rbindlist(out$best_tune, idcol="source")
  out$caret_train <- full_data_model
  out$seed <- seed
  
  return(out)
}

