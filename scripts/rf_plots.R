library(dplyr)
library(data.table)
library(tidyr)
library(precrec)
library(stringr)
library(tibble)

library(ggplot2)
library(ggpubr)
library(lemon)
library(drlib)
library(ggrepel)

###############################################################
# some utility functions
###############################################################

RANK_ORDER <- c("Class", "Order", "Family", "Genus", "Species")

get_main_text_rows <- function(df) {
  # returns the experiment settings used in the main text
  df %>%
    filter(
      transform=="relative_abundance",
      min_taxa_prev==20,
      k==5,
      source=="observed"
    )
}

get_pvalue_label <- function(pvals) {
  out <- rep(NA, length(pvals))
  out[ pvals < 0.1 ] <- "*"
  out[ pvals < 0.05 ] <- "**"
  return(out)
}

plot_theme <- function() {
  theme(
    strip.text=element_text(colour="black"),
    strip.background=element_rect(fill=NA, colour=NA)
  )
}

###############################################################
###############################################################

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
# make all panels for Figure 1
###############################################################

combined_taxonomy <- lapply(
  phylos,
  function(x) x@tax_table@.Data %>% 
    as_tibble()
) %>%
  bind_rows()

save_path <- "../results/main_text"

# AUCs and p-values
aucs_w_pvalues <- arrow::read_feather(
  file.path(save_path, "performance_aucs_pvalues.feather")
) %>%
  rename(type=curvetypes) %>%
  group_by(phenotype, transform, k, type) %>%
  mutate(p_value=p.adjust(p_value, method="fdr")) %>%
  ungroup() %>%
  mutate(kingdom=recode(kingdom, all="both") %>%
           factor(levels=c("bacterial", "fungal", "both")))

# AU ROC panels
PVAL_AST_YPOS <- 1.0
tick_font_size <- 10
axis_title_font_size <- 12

auc_plot_titles <- c(
  Disease="BX vs CF",
  Group="FB vs NAFD (CF only)",
  Exacing="CFPE vs No CFPE (CF only)"
)

auroc_panels <- aucs_w_pvalues %>%
  filter(!phenotype %in% c("Class.fungal", "Genus.bacterial")) %>%
  get_main_text_rows() %>%
  filter(type=="ROC") %>%
  mutate(pvalue_label=get_pvalue_label(p_value)) %>%
  group_by(phenotype) %>%
  do(
    plot=ggplot(
      data=.,
      aes(x=factor(glom_rank, levels=RANK_ORDER), y=auc, fill=kingdom)
      ) +
      geom_col(position="dodge") +
      geom_text(aes(label=pvalue_label, y=upper_bound+0.02), 
                position=position_dodge(0.9), size=5) +
      geom_hline(yintercept=0.5,
                 linetype="dashed", colour="red", size=1) +
      geom_errorbar(aes(ymin=lower_bound, ymax=upper_bound),
                    position=position_dodge(0.9)) +
      scale_y_continuous(limits=c(0,1.02), breaks=seq(0,1,by=0.25)) +
      ylab("AU-ROC") +
      xlab("Rank") +
      labs(fill="Kingdom") +
      theme_linedraw() +
      plot_theme() +
      theme(
        panel.grid.major=element_line(size=0.01, colour="grey"),
        panel.grid.minor=element_line(size=0.01, colour="grey"),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.text=element_text(size=tick_font_size),
        axis.title=element_text(size=axis_title_font_size),
        legend.text=element_text(size=tick_font_size),
        legend.title=element_text(size=axis_title_font_size)
      ) +
      ggtitle(auc_plot_titles[[ .$phenotype[[1]] ]])
  )
auroc_panels <- setNames(auroc_panels$plot, auroc_panels$phenotype)

# AU PRC panels
prc_baseline_values <- lapply(sample_groups[1:3], function(x) x %>% droplevels %>% as.integer()-1)
prc_baseline_values <- sapply(prc_baseline_values, function(x) sum(x)/length(x)) 

auprc_panels <- aucs_w_pvalues %>%
  filter(!phenotype %in% c("Class.fungal", "Genus.bacterial")) %>%
  get_main_text_rows() %>%
  filter(type=="PRC") %>%
  mutate(pvalue_label=get_pvalue_label(p_value)) %>%
  group_by(phenotype) %>%
  do(
    plot=ggplot(
      data=.,
      aes(x=factor(glom_rank, levels=RANK_ORDER), y=auc, fill=kingdom)
    ) +
      geom_col(position="dodge") +
      geom_text(aes(label=pvalue_label, y=upper_bound+0.02), 
                position=position_dodge(0.9), size=5) +
      geom_hline(yintercept=prc_baseline_values[[ .$phenotype[[1]] ]],
                 linetype="dashed", colour="red", size=1) +
      geom_errorbar(aes(ymin=lower_bound, ymax=upper_bound),
                    position=position_dodge(0.9)) +
      scale_y_continuous(limits=c(0,1.02), breaks=seq(0,1,by=0.25)) +
      ylab("AU-PRC") +
      xlab("Rank") +
      labs(fill="Kingdom") +
      theme_linedraw() +
      plot_theme() +
      theme(
        panel.grid.major.y=element_line(size=0.01, colour="grey"),
        panel.grid.minor.y=element_line(size=0.01, colour="grey"),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.text=element_text(size=tick_font_size),
        axis.title=element_text(size=axis_title_font_size),
        legend.text=element_text(size=tick_font_size),
        legend.title=element_text(size=axis_title_font_size)
      ) +
      ggtitle(auc_plot_titles[[ .$phenotype[[1]] ]])
  )
auprc_panels <- setNames(auprc_panels$plot, auprc_panels$phenotype)

auc_fig <- ggarrange(
  auroc_panels$Disease, auprc_panels$Disease,
  auroc_panels$Group, auprc_panels$Group,
  auroc_panels$Exacing, auprc_panels$Exacing,
  ncol=2, nrow=3,
  common.legend=TRUE,
  legend="right",
  labels="auto"
)
ggsave("../plots/rf_auc_fig.png", dpi=450, plot=auc_fig)


###############################################################
###############################################################

###############################################################
# make all panels for Figure 1
# Table S1:AUCs under different transformations, values of k
# and minimum prevalence threshold
###############################################################

library(writexl)

auc_supp_tab_data <- aucs_w_pvalues %>%
  select(-c(n, spec_name, n_perms, error, lower_bound, upper_bound)) %>%
  filter(!phenotype %in% c("Class.fungal", "Genus.bacterial"), source=="observed") %>%
  janitor::remove_constant() %>%
  mutate(auc_pval=sprintf("%.2f (%.2f)", auc, p_value)) %>%
  select(-c(auc, p_value)) %>%
  mutate(column_title=sprintf("minimum prevalence=%d%%", min_taxa_prev)) %>%
  select(-c(min_taxa_prev)) %>%
  pivot_wider(names_from=column_title, values_from=auc_pval)

tbl_list <- expand_grid(
  phenotype=unique(auc_supp_tab_data$phenotype),
  kingdom=unique(auc_supp_tab_data$kingdom)
)
tbl_list$dat <- apply(
  tbl_list,
  1,
  function(x) auc_supp_tab_data %>%
    filter(phenotype==x[[1]], kingdom==x[[2]]) %>%
    select(-c(phenotype, kingdom))
  )
  
setNames(
  tbl_list$dat,
  apply(tbl_list, 1, function(x) paste(x[[1]], "-", x[[2]]))
) %>%
  write_xlsx(path="../csv/supplementary_table_1.xlsx")

###############################################################
###############################################################

###############################################################
# permutation p-values for comparuing the AUCs
# of models trained using each (or both) kingdom
###############################################################

aucs_with_perms <- arrow::read_feather(
  file.path(save_path, "performance_aucs_with_perms.feather")
)

auc_diff_testdata <- aucs_with_perms %>%
  select(-c(error, lower_bound, upper_bound, n)) %>%
  pivot_wider(names_from=source, values_from=auc) %>%
  group_by(phenotype, transform, glom_rank, k, curvetypes, observed, kingdom, min_taxa_prev) %>%
  do(permuted_aucs=.[,11:ncol(.)] %>% slice(1) %>% as.numeric())

perm_tbl <- auc_diff_testdata %>%
  select(-observed) %>%
  group_by(phenotype, transform, glom_rank, k, curvetypes) %>%
  tidybayes::gather_pairs(key=kingdom, value=permuted_aucs, triangle="both only")
perm_tbl$perm_auc_diff <- mapply(function(x,y) x-y, perm_tbl$.x, perm_tbl$.y, SIMPLIFY=FALSE)

obs_tbl <- auc_diff_testdata %>%
  select(-permuted_aucs) %>%
  group_by(phenotype, transform, glom_rank, k, curvetypes) %>%
  tidybayes::gather_pairs(key=kingdom, value=observed, triangle="both only") %>%
  mutate(obs_auc_diff=.x-.y)

perm_pvalues <- perm_tbl %>%
  select(-c(.y, .x)) %>%
  inner_join(
    obs_tbl 
  ) 

calc_pvalue <- function(obsval, perm_val_list, alternative=c("less", "greater")) {
  # compute a permutation p-value given observed AUC differences and permuted values
  alternative <- match.arg(alternative)
  if(alternative=="greater") {
    n_greater <- sum(obsval > perm_val_list, na.rm=TRUE)
    n_valid <- sum(!is.na(perm_val_list))
    return((1+n_greater)/(n_valid+1))
  } else {
    n_less <- sum(obsval < perm_val_list, na.rm=TRUE)
    n_valid <- sum(!is.na(perm_val_list))
    return((1+n_less)/(n_valid+1))
  }
}

perm_pvalues$pvalue_less <- mapply(calc_pvalue, perm_pvalues$obs_auc_diff, perm_pvalues$perm_auc_diff,
                                   alternative="less")
perm_pvalues$pvalue_greater <- mapply(calc_pvalue, perm_pvalues$obs_auc_diff, perm_pvalues$perm_auc_diff,
                                   alternative="greater")
perm_pvalues$n_valid <- sapply(perm_pvalues$perm_auc_diff, function(x) sum(!is.na(x)))

perm_pvalues <- perm_pvalues %>%
  pivot_longer(c(pvalue_less, pvalue_greater), names_to="alternative", values_to="pvalue") %>%
  mutate(alternative=str_remove(alternative, "pvalue_"))

# permutation p-values for adding second kingdom to RF model
# corrected using FDR
perm_pvalues %>%
  ungroup() %>%
  filter(transform=="relative_abundance", k==5, min_taxa_prev==20, alternative=="greater") %>%
  select(phenotype, glom_rank, curvetypes, .row, .col, pvalue) %>%
  mutate(pvalue=p.adjust(pvalue, method="fdr")) %>%
  pivot_wider(names_from=glom_rank, values_from=pvalue) %>%
  filter(.row=="all", phenotype %in% c("Disease", "Group")) %>%
  arrange(phenotype)

# p-values for adding bacteria to model
perm_pvalues %>%
  ungroup() %>%
  filter(transform=="relative_abundance", k==5, min_taxa_prev==20, alternative=="less") %>%
  select(phenotype, glom_rank, curvetypes, .row, .col, pvalue) %>%
  mutate(pvalue=p.adjust(pvalue, method="fdr")) %>%
  pivot_wider(names_from=glom_rank, values_from=pvalue) %>%
  filter(.row=="bacterial", phenotype %in% c("Disease", "Group")) %>%
  arrange(phenotype)

###############################################################
###############################################################

###############################################################
# RF variable importance barplots (Figure 2 a,c,e,g)
###############################################################

taxa_renamer <- setNames(combined_taxonomy$Genus, combined_taxonomy$OTUID)
kingdom_getter <- setNames(combined_taxonomy$Domain, combined_taxonomy$Genus)

varimp_values <- arrow::read_feather(
  file.path(save_path, "combined_var_imps_observed.feather")
) %>%
  as_tibble() %>%
  group_by(phenotype, kingdom, transform, spec_name) %>%
  mutate(pvalue=p.adjust(pvalue, method="fdr")) %>%
  ungroup()

varimp_barplot_data <- varimp_values %>%
  get_main_text_rows() %>%
  group_by(phenotype, kingdom, spec_name) %>%
  filter(!grepl("conf", variable)) %>%
  slice_max(importance, n=5) %>%
  ungroup() %>%
  mutate(
    new_variable=taxa_renamer[ variable ],
    taxa_kingdom=factor(
      kingdom_getter[ new_variable ],
      levels=c("k__Bacteria", "k__Fungi"),
      labels=c("Bacteria", "Fungi")
    ),
    new_variable=str_remove(new_variable, "g__|f__"),
    pvalue_label=get_pvalue_label(pvalue),
    imp_method=str_extract(spec_name, "impurity_corrected|impurity|permutation")
  )

pretty_facet_titles <- c(
  "Disease--all"="CF vs BX (both kingdoms)",
  "Disease--bacterial"="CF vs BX (bacteria)",
  "Disease--fungal"="CF vs BX (fungi)",
  "Group--fungal"="FB vs NAFD (CF only, fungi)"
)

varimp_barplot <- varimp_barplot_data %>%
  mutate(
    new_variable=if_else(grepl("___conf", variable), str_remove(variable, "___conf"), new_variable),
    imp_method=factor(
      imp_method,
      levels=c("permutation", "impurity", "impurity_corrected"),
      labels=c("MDA", "MDG", "Corrected MDG")
    ),
    model_name=pretty_facet_titles[ sprintf("%s--%s", phenotype, kingdom) ]
  ) %>%
  group_by(phenotype, kingdom) %>%
  do(
    plot=ggplot(data=.,
                aes(x=reorder_within(new_variable, -importance, imp_method),
                    y=importance,
                    fill=taxa_kingdom)) +
      geom_col() +
      geom_text(aes(label=pvalue_label), vjust=0.15,
                size=5) +
      facet_wrap(~imp_method, scales="free") +
      theme_linedraw() +
      plot_theme() +
      scale_x_reordered() +
      ylab("Importance") +
      xlab("Taxa") +
      theme(
        axis.text.x=element_text(angle=30, vjust=1, hjust=1),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(size=0.01, colour="grey"),
        panel.grid.minor.y=element_line(size=0.01, colour="grey"),
        strip.text=element_text(size=12),
        axis.title.x=element_blank()
      ) +
      labs(fill="Kingdom") +
      scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
      scale_fill_manual(drop=FALSE, values=scales::hue_pal()(3)[1:2])
      
  )

ggarrange(
  plotlist=varimp_barplot$plot,
  ncol=1, nrow=4,
  labels=c("auto"),
  common.legend=TRUE
)
ggsave("../plots/rf_varimp_barplotsw_confounding.pdf", height=9, width=6)

##########################################################
##########################################################

##########################################################
# partial dependence plots (Figure 2 b,d,f,h)
##########################################################

partial_dep_data <- lapply(
  list.files(file.path(save_path, "pdp_data"), full.names=TRUE),
  function(x) fread(x) %>% mutate(x=as.numeric(x))
) %>%
  bind_rows() %>%
  as_tibble() %>%
  filter(!is.na(x))

y_labels <- c(
  Disease=sprintf('RF prediction\n(\u27F5 CF)\t(BX \u27F6)'),
  Group=sprintf(' RF prediction\n(\u27F5 NAFD)\t(FB \u27F6)')
)

pretty_facet_titles <- c(
  "Disease--all"="CF vs BX (both kingdoms)",
  "Disease--bacterial"="CF vs BX (bacteria)",
  "Disease--fungal"="CF vs BX (fungi)",
  "Group--fungal"="FB vs NAFD (CF only, fungi)"
)

pdp_panel_plotdata <- partial_dep_data %>%
  inner_join(varimp_barplot_data) %>%
  filter(!grepl("___conf", variable)) %>%
  filter(transform=="relative_abundance") %>%
  mutate(model_name=pretty_facet_titles[ sprintf("%s--%s", phenotype, kingdom) ]) 

colour_palette <- setNames(
  scales::hue_pal()(length(unique(pdp_panel_plotdata$new_variable))),
  unique(pdp_panel_plotdata$new_variable)
)

pdp_panels <- pdp_panel_plotdata %>%
  group_by(phenotype, kingdom) %>%
  do(
    plot=ggplot(data=., aes(x=x, y=yhat, colour=new_variable)) +
      geom_line(aes(group=new_variable), size=0.75) +
      theme_linedraw() +
      plot_theme() +
      geom_text_repel(data=. %>% group_by(new_variable) %>% slice_sample(n=1),
                      aes(label=new_variable), show.legend=FALSE) +
      facet_wrap(~model_name) +
      theme(legend.position="none",
            axis.text=element_text(size=10),
            axis.title=element_text(size=12),
            strip.text=element_text(size=12),
            panel.grid.major=element_line(size=0.01, colour="grey"),
            panel.grid.minor=element_line(size=0.01, colour="grey")) +
      scale_y_continuous(limits=c(-3.5, 0.5)) +
      scale_colour_manual(values=colour_palette) +
      ylab(y_labels[[ .$phenotype[[1]] ]]) +
      xlab("Relative abundance")
  )

# combine RF importances and PDP panels into a single plot (Figure 2)
all_panels <- list(
  lapply(varimp_barplot$plot, function(x) x + theme(legend.position="none")),
  pdp_panels$plot) %>%
  unlist(recursive = FALSE)

fig <- cowplot::plot_grid(
  plotlist=all_panels,
  ncol=2, nrow=4,
  align="hv", axis="lrtb",
  rel_widths=c(1.3, 1),
  byrow=FALSE,
  labels="auto"
)

ggsave("../plots/rf_varimp_fig2.svg", plot=fig, height=9, width=9)

##########################################################
##########################################################

##########################################################
# Variable importance values for different transformations
##########################################################

library(GGally)

varimp_pairs_plots <- varimp_values %>%
  filter(!grepl("__conf", variable)) %>%
  filter(source=="observed") %>%
  janitor::remove_constant() %>%
  mutate(imp_method=str_extract(spec_name, "impurity_corrected|impurity|permutation")) %>%
  select(-c(pvalue, spec_name)) %>%
  pivot_wider(names_from=imp_method, values_from=importance) %>%
  group_by(phenotype, kingdom, min_taxa_prev) %>%
  do(
    plot=GGally::ggpairs(data=.,
                         columns=c("impurity_corrected", "impurity", "permutation"),
                         aes(colour=transform),
                         upper=list(continuous=wrap("cor", method="spearman", size=3)),
                         lower=list(continuous=wrap("points", size=1.5)),
                         switch="both") +
      theme_linedraw() +
      plot_theme() +
      theme(strip.placement="outside",
            strip.text=element_text(size=11),
            plot.margin=margin(10, 10, 10, 10, "pt"))
  )

for(phenotype_name in unique(varimp_pairs_plots$phenotype)) {
  for(kingdom_name in unique(varimp_pairs_plots$kingdom)) {
    plot_tbl <- varimp_pairs_plots %>%
      filter(phenotype==phenotype_name, kingdom==kingdom_name)
    if(nrow(plot_tbl)==0) next
    
    fig <- ggarrange(
      plotlist=lapply(plot_tbl$plot, ggmatrix_gtable),
      labels="auto",
      label.y=1.01
    ) + bgcolor("White")
    
    ggsave(
      sprintf("../plots/varimp_ggpairs_%s_%s3.png", phenotype_name, kingdom_name),
      plot=fig,
      width=10, height=10
    )
  }
}

##########################################################
##########################################################
