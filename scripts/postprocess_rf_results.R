library(pbmcapply)
library(data.table)
library(dplyr)
library(multidplyr)
library(tidyr)
library(stringr)
library(precrec)


if(interactive()) {
  N_CORES <- 2
  save_path <- c("../results/DEBUG")
} else {
  N_CORES <- Sys.getenv("N_CORES") %>% as.integer()
  save_path <- Sys.getenv("SAVE_PATH")
}

message("save_path: ", save_path)
message("N_CORES: ", N_CORES)


get_perf_curves <- function(probs, labs, fold_ids, alpha=0.95, x_bins=20) {
  
  out <- tryCatch({
    
    calc_auc_cis <- !all(is.na(fold_ids))
    
    if(!calc_auc_cis) {
      model_eval_data <- precrec::mmdata(
        scores=probs, labels=labs
      )
    } else {
      model_eval_data <- precrec::mmdata(
        nfold_df=data.frame(score=probs, lab=labs, fold=fold_ids),
        score_col="score",
        lab_col="lab",
        fold_col="fold")
    }
    
    basic_eval_plotdata <- model_eval_data %>%
      precrec::evalmod(mode="basic") %>%
      tidyr::as_tibble() %>%
      dplyr::select(-modname)
    
    curve_eval <- model_eval_data %>%
      precrec::evalmod(cb_alpha=1-alpha,
                       x_bins=x_bins)
    
    if(calc_auc_cis) {
      auc_df <- curve_eval %>%
        precrec::auc_ci(alpha=1-alpha) %>%
        dplyr::select(-modnames) %>%
        dplyr::rename(auc=mean)
      
      curve_data <- curve_eval %>%
        tidyr::as_tibble() %>%
        dplyr::select(-modname)
    } else {
      auc_df <- curve_eval %>%
        precrec::auc() %>%
        dplyr::select(-c(dsids, modnames)) %>%
        dplyr::rename(auc=aucs)
      
      curve_data <- tibble::tibble()
    }
    
    all_curve_data <- basic_eval_plotdata %>%
      dplyr::bind_rows(curve_data)
    
    return(list(curve_data=all_curve_data, auc_df=auc_df))
  },
  error=function(cond) {
    return(list(curve_data=tibble::tibble(), auc_df=tibble::tibble()))
  }
  )
  
  return(out)
}


# concatenate all the relevant results
for(dir_name in c("validation_fold_preds", "var_imps")) {
  cat(dir_name, "\n") 
  
  csv_names <- lapply(
    save_path,
    function(x) list.files(file.path(x, dir_name), full.names=TRUE)
  ) %>% unlist(recursive=FALSE)
  stopifnot(length(unique(csv_names))==length(unique(basename(csv_names))))
  names(csv_names) <- csv_names %>% basename()
  cat(sprintf("Found %s files\n", length(csv_names) %>% scales::comma()))
  
  if(length(csv_names)>0) {
  
    # concatenate all the validation fold preds from all runs and save them as a
    # single dataframe
    tmp <- pbmclapply(
      csv_names, 
      fread,
      mc.cores=N_CORES
    ) %>% rbindlist(use.names=TRUE, fill=TRUE)
    
    cat(sprintf("Saving to %s\n", file.path(save_path, sprintf("combined_%s.feather", dir_name))))
    arrow::write_feather(tmp, file.path(save_path, sprintf("combined_%s.feather", dir_name)))
    
  } else {
    # nothing to do here for variable importance
    warning("No csv files found", call.=FALSE)
    next
  }

  # compute the performance metrics
  if(dir_name=="validation_fold_preds") {
    cat("Computing predictive performance metrics\n")
    
    cl <- new_cluster(N_CORES) %>%
      cluster_library(c("dplyr")) %>%
      cluster_copy("get_perf_curves")

    # ROC and PR curve data (not just AUCs)
    cat("Computing ROC and PR curves...")
    perf_result <- tmp %>%
      mutate(outer_fold=if_else(source=="observed", outer_fold, NA_integer_)) %>%
      group_by_at(
        setdiff(colnames(.), c("outer_fold", "predictions", "labels", "sample_id"))
      ) %>%
      partition(cl) %>%
      summarise(
        curve_result=list(get_perf_curves(
          probs=predictions, labs=labels, fold_ids=outer_fold))
        ) %>%
      collect()
      # unnest(curve_data)
    
    perf_result$curve_data = lapply(perf_result$curve_result, function(x) x[[1]])
    perf_result$auc_df = lapply(perf_result$curve_result, function(x) x[[2]])
    
    perf_curves <- perf_result %>%
      unnest(curve_data) %>%
      select(-c(curve_result, auc_df)) %>%
      filter(source=="observed")
    
    # save the curve data
    perf_curves %>%
      arrow::write_feather(file.path(save_path, "performance_metric_results.feather"))
    
    cat("done\n")
    
    # compute statistical significance of AUCs
    perf_aucs <- perf_result %>%
      unnest(auc_df) %>%
      select(-c(curve_data, curve_result))
    perf_aucs %>%
      arrow::write_feather(file.path(save_path, "performance_aucs_with_perms.feather"))
    
    # p-values of AUCs from label randomisation test
    cat("Computing permutation p-values...")
    auc_pvalues <- perf_aucs %>%
      select(-c(error, lower_bound, upper_bound, n)) %>%
      pivot_wider(names_from=source, values_from=auc)  %>%
      pivot_longer(contains("permutation_"), names_to="resample_name", values_to="resample_value") %>%
      group_by_at(
        setdiff(colnames(.), c("resample_name", "resample_value"))
        ) %>%
      summarise(
        n_greater=sum(observed<resample_value, na.rm=TRUE),
        n_reps=sum(!is.na(resample_value)),
        .groups="keep"
        ) %>%
      ungroup() %>%
      mutate(p_value=(1+n_greater)/(1+n_reps)) %>%
      select(-c(n_greater, n_reps)) %>%
      rename(auc=observed)
    cat("done\n")
    
    cat("Saving...")
    perf_aucs <- perf_aucs %>% 
      inner_join(auc_pvalues)

    perf_aucs %>%
      arrow::write_feather(file.path(save_path, "performance_aucs_pvalues.feather"))
    cat("done\n")
    
  }
}

cat("Script finished successfully\n")

