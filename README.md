# Cross-kingdom analysis of microbial communities in Cystic Fibrosis and Bronchiectasis

The repo contains the code to generate the results in [Cuthberton and Ish-Horowicz, (2022)](https://www.biorxiv.org/content/10.1101/2022.01.11.475678v2).

# Generating results

All code assumes that `run` is the current working directory.

## Preprocessing

The processing of both the bacterial and fungal QIIME files into phyloseqs can be found in [markdown/FAME_preprocessing.Rmd](markdown/FAME_preprocessing.Rmd). Note that this script is not runnable as the patient metadata files are not included to ensure data protection. Phyloseq objects that have been stripped of their full metadata (plus the metadata columns required for the analysis) are available in the `data` directory.

## Random forest experiments (Figures 1 and 2)

### Running the experiments

The random forest two-sample testing script should be called via `run/run_all_rf_experiments.sh`. There are two types of experiment - the type depends on the environment variable `VARIMP_RUN` (set within the script):

* `VARIMP_RUN=FALSE`: estimates the area under curves via nested cross-valdiation (Figure 1).

* `VARIMP_RUN=TRUE`: fits the model on the full datasets and computes variable importance scores (Figure 2).

The script calls [scripts/run_rf_experiments.R](scripts/run_rf_experiments.R).

The label randomisation procedure used to assess the statistical significance of the area under curves requires fitting a large number of models. The validation (held-out) predictions for each  model are saved as separate csv files in the `results/PBS_JOBID` directory, where `PBS_JOBID` is an environment variable set in `run/run_all_rf_experiments.sh`. Once all these results have been completed, the AUCs and corresponding p-values are computed by `run/compute_rf_aucs.sh`, where you should set the `SAVE_PATH` environment variable to `results/PBS_JOBID` so it knows where to find the csv files.

This generates four files:

* `combined_validation_fold_preds.feather`: the validation set predictions from all models
* `performance_aucs_pvalues.feather`: the AUCs on the observed and the corresponding p-values from the label randomisation test
* `performance_aucs_with_perms.feather`: the AUCs with those from the label randomisation test
* `performance_metric_results.feather`: the coordinates of the PR and ROC curves

The paper contains many different configurations under which the random forests are trained. If `VARIMP_RUN=FALSE` the configurations are

* three host phenotypes
* three sets of covariates (bacteria, fungi and both)
* five transformations
* five taxonomic ranks for agglomeration
* two values of $k$ (5,10) in the $k$-fold cross-validation 
* three prevalence thresholds for removing rare taxa (0%, 10%, 20%)

This is 5,400 different configurations. The `PBS_ARRAY_INDEX` argument controls which configuration is run so is easy to run them all using an array job.

If `VARIMP_RUN=TRUE` then there are only 180 configurations as we only agglomerate to Genus and only include models that are able to discriminate between the classes. The value of $k$ is also ignored since cross-validation is not required to compute the variable importance scores.

### Generating the plots

The code used to generate Figures 1 and 2 is in [scripts/rf_plots.R](scripts/rf_plots.R). All the required plot data can be found in [results/main_text](results/main_text).

## Unsupervised analysis (Figure 3)

### Network analysis using SparCC (Figures 3a and 3b)

Run the [scripts/network_analysis.R](scripts/network_analysis.R) script to generate the main text heamap of SparCC correlations and the boxplots comparing inter- and intra-kingdom correlations.

### MOFA (Figure 3c)

Run the [scripts/mofa_analysis.R](scripts/mofa_analysis.R) script to generate the plot showing the variance explained by the MOFA factors at each of the taxonomic ranks.

## Dirichlet Multinomial Mixture clustering (Figure 4)

Run the [scripts/dmn_clustering.R](scripts/dmn_clustering.R) script to run the sample-wise clustering analysis and generate all panels of Figure 4.

## Ecological analyses



# Relevant citations

Cuthbertson, L., Ish-Horowicz, J., Felton, I., James, P., Turek, E., Cox, M. J., ... & Cookson, W. O. (2022). Cross-kingdom analysis of microbial communities in Cystic Fibrosis and Bronchiectasis. bioRxiv.

Cuthbertson, L., Felton, I., James, P., Cox, M. J., Bilton, D., Schelenz, S., ... & Moffatt, M. F. (2021). The fungal airway microbiome in cystic fibrosis and non-cystic fibrosis bronchiectasis. Journal of Cystic Fibrosis, 20(2), 295-302.