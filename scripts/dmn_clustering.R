library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(DirichletMultinomial)
library(ggplot2)
library(ggpubr)
library(ggforce)

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
# prepare OTU tables
###############################################################

combined_taxonomy <- lapply(
  phylos,
  function(x) x@tax_table@.Data %>% 
    as_tibble()
  ) %>%
  bind_rows()

TAXA_ORDER <- c("Class", "Order", "Family", "Genus", "Species")

phylos_glom <- lapply(
  setNames(TAXA_ORDER, TAXA_ORDER),
  function(rank_name) lapply(phylos, function(x) speedyseq::tax_glom(x, rank_name))
) 

X <- lapply(
  phylos_glom,
  function(x) lapply(x, function(xx) xx %>% otu_table() %>% as.data.frame() %>% t())
) %>% unlist(recursive=FALSE)
X <- lapply(X, function(x) x[ rownames(y), ])
lapply(X, dim)

###############################################################
###############################################################

###############################################################
# run DMM analysis
###############################################################

get_gof <- function(fit) {
  return(c("laplace"=laplace(fit), "AIC"=AIC(fit), "BIC"=BIC(fit)))
}

max_k <- 5 # maximum number of clusters
SEED <- 123412345
set.seed(SEED)

dmn_fits <- pbmcapply::pbmclapply(
  X,
  function(x) lapply(
    1:max_k,
    function(k) dmn(x, k=k, seed=SEED)
  ),
  mc.cores=10
)

# check goodness of fits using different numbers of clusters
dmn_fit_vals <- lapply(dmn_fits,
                       function(x) lapply(x,
                                          function(xx) goodnessOfFit(xx))
                       %>% bind_rows() %>%
                         rownames_to_column("k")
) %>%
  bind_rows(.id="kingdom") %>%
  separate(kingdom, c("rank", "kingdom"), sep="\\.") %>%
  pivot_longer(-c(kingdom, k, rank)) %>%
  filter(name %in% c("AIC", "BIC", "Laplace")) %>%
  mutate(k=as.integer(k),
         kingdom=recode(kingdom, bacterial="Bacteria", fungal="Fungi"))

dmn_fit_vals %>%
  filter(k<=5) %>%
  ggplot(aes(x=k, y=value, colour=rank,
             group=interaction(name, rank))) +
  geom_point() +
  geom_line() +
  lemon::facet_rep_grid(rows=vars(kingdom), cols=vars(name),
                        scales="free_y") +
  theme_classic(base_size=16) +
  xlab("Number of clusters") + 
  ylab("Goodness of fit metric value") +
  geom_point(
    data=dmn_fit_vals %>%
      group_by(rank, kingdom, name) %>%
      summarise(k=k[ which.min(value) ], value=value[ which.min(value) ]),
    shape=21, size=5, stroke=1
  ) +
  labs(colour="Rank") +
  theme(panel.spacing=unit(0, "lines"))

#
# check cluster robustness using consensus clustering
# - how often are taxa clustered together in a bootstrap replicate (proprtionality)?
library(ClusBoot)
library(ComplexHeatmap)

do_dmn_clust <- function(x, k, ...) {
  x %>% dmn(k=k) %>% mixture(assign=TRUE)
}

clustboot_fungal <- clusboot(X$Class.fungal, B=100, clustering.func=purrr::partial(do_dmn_clust, k=2))
clustboot_bacterial <- clusboot(X$Genus.bacterial, B=100, clustering.func=purrr::partial(do_dmn_clust, k=2))

# create a heatmap showing the proportionality
propr_heatmap <- function(res, title, cluster_colours) {
  Heatmap(
    res$proportions,
    name="proportion",
    col=circlize::colorRamp2(c(0,1), c("white", "darkgreen")),
    bottom_annotation=columnAnnotation(cluster=as.factor(res$clustering),
                                       col=cluster_colours,
                                       show_legend=FALSE),
    right_annotation=rowAnnotation(cluster=as.factor(res$clustering),
                                   col=cluster_colours,
                                   show_legend=FALSE,
                                   show_annotation_name=FALSE),
    cluster_rows=FALSE, cluster_columns=FALSE,
    show_row_names=FALSE, show_column_names=FALSE,
    row_names_gp=gpar(fontsize=7),
    column_names_gp=gpar(fontsize=7),
    column_title=title,
    column_title_gp=gpar(fontsize=20)
  )
}

hm_fungal <- clustboot_fungal %>% 
  propr_heatmap("Fungal (Class) DMN clusters",
                cluster_colours=list(cluster=c("1"="blue", "2"="red")))
hm_bacterial <- clustboot_bacterial %>% 
  propr_heatmap("Bacterial (Genus) DMN clusters",
                cluster_colours=list(cluster=c("1"="purple", "2"="orange")))

png("../plots/submission/figureS7.png", width=16, height=8, units="in", res=600)
draw(hm_fungal + hm_bacterial, ht_gap=unit(1, "cm"), auto_adjust=FALSE, merge_legend=TRUE)
dev.off()




# get cluster assignemtns from the best fits
best_fits <- list(Class.fungal=dmn_fits$Class.fungal[[2]],
                  Genus.bacterial=dmn_fits$Genus.bacterial[[2]])
dmm_clusters <- lapply(
  best_fits,
  function(x) mixture(x, assign=TRUE) %>% as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    rename("."="cluster")
) %>%
  bind_rows(.id="kingdom") %>%
  pivot_wider(names_from=kingdom, values_from=cluster)


# adjusted Rand index between cluster assignments and
# host phenotypes
dmm_clusters %>%
  left_join(
    y %>% rownames_to_column("sample_id"),
    by="sample_id"
  ) %>%
  select(-sample_id) %>%
  corrr::colpair_map(mclust::adjustedRandIndex)

#
# make PCoA plots
library(vegan)
library(ape)

get_biplot_data <- function(x, plot.axes=c(1,2)) {
  pr.coo <- x$vectors
  diag.dir <- diag(c(1,1))
  pr.coo[,plot.axes] <- pr.coo[,plot.axes] %*% diag.dir
  return(pr.coo[,plot.axes] %>%
           as.data.frame() %>%
           rownames_to_column("sample_id"))
}

bc_dists <- lapply(X[ names(best_fits) ], function(x) vegdist(x, method="bray"))
pcoa_res <- lapply(bc_dists, ape::pcoa)

biplot_data <- lapply(pcoa_res, get_biplot_data)

axis_labels <- list(
  "Fungi (Class)"=c("PCoA axis 1 (32%)", "PCoA axis 2 (21%)"),
  "Bacteria (Genus)"=c("PCoA axis 1 (38%)", "PCoA axis 2 (11%)")
)

biplots <- biplot_data %>%
  bind_rows(.id="kingdom") %>%
  left_join(dmm_clusters %>%
              pivot_longer(-sample_id, names_to="kingdom", values_to="cluster"),
            by=c("kingdom", "sample_id")) %>%
  mutate(cluster=factor(paste0(kingdom, cluster),
                        levels=c("Genus.bacterial1", "Genus.bacterial2", "Class.fungal1", "Class.fungal2"),
                        #labels=c("Bacteria1", "Bacteria2", "Fungi1", "Fungi2")
                        labels=c("No Pseudomonas\ndomination", "Pseudomonas\ndomination",
                                 "Saccharomycetes\ndomination", "No Saccharomycetes\ndomination")),
         kingdom=recode(kingdom, Class.fungal="Fungi (Class)", Genus.bacterial="Bacteria (Genus)")) %>%
  group_by(kingdom) %>%
  do(
    plot=ggplot(data=., aes(x=Axis.1,
                            y=Axis.2,
                            colour=as.factor(cluster))) +
      geom_point() +
      # scale_colour_discrete() +
      xlab(axis_labels[[.$kingdom[[1]]]][[1]]) +
      ylab(axis_labels[[.$kingdom[[1]]]][[2]]) +
      theme_linedraw() +
      theme_classic() +
      stat_ellipse() +
      # geom_mark_ellipse(aes(label=as.factor(cluster))) +
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=12),
            strip.text=element_text(size=12),
            panel.grid.major=element_line(size=0.01, colour="grey"),
            panel.grid.minor=element_line(size=0.01, colour="grey"),
            aspect.ratio=1) +
      labs(colour="Cluster") +
      coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))
  )

#
# stacked barplots with samples ordered by cluster label
taxa_renamer <- list(
  Class.fungal=setNames(combined_taxonomy$Class, combined_taxonomy$OTUID),
  Genus.bacterial=setNames(combined_taxonomy$Genus, combined_taxonomy$OTUID)
) 

displayed_taxa <- lapply(
  setNames(names(best_fits), names(best_fits)),
  function(x) tibble(taxa=taxa_renamer[[x]][ colnames(X[[x]]) ],
                     reads=colSums(X[[x]]))) %>%
  bind_rows(.id="kingdom") %>%
  group_by(kingdom) %>%
  slice_max(order_by=reads, n=5) %>%
  ungroup()

make_colour_map <- function(x) {
  unique_x <- unique(x)
  cmap <- scales::hue_pal()(length(unique_x)-1)
  names(cmap) <- unique_x[ unique_x!="Other" ]
  cmap[["Other"]] <- "grey50"
  return(cmap)
}

# create plots (one for each kingdom)
stacked_barplots <- list()

for(plot_name in c("Class.fungal", "Genus.bacterial")) {
  
  plot_data <- X[[ plot_name ]] %>%
    relative_abundance() %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(-sample_id) %>%
    left_join(dmm_clusters, by="sample_id") %>%
    mutate(name=taxa_renamer[[plot_name]][ name ]) %>%
    ungroup() %>%
    mutate(
      name=if_else(name %in% displayed_taxa$taxa, name, "Other"),
      Class.fungal=recode(Class.fungal, "1"="Saccharomycetes domination", "2"="No Saccharomycetes domination"),
      Genus.bacterial=recode(Genus.bacterial, "1"="No Pseudomonas domination", "2"="Pseudomonas domination")
    )
  
  stacked_barplots[[plot_name]] <- plot_data %>%
    ggplot(aes(x=sample_id, y=value, fill=name)) +
    geom_col(colour="black", size=0.25) +
    ggforce::facet_row(as.formula(paste0("~", plot_name)), scales="free_x", space="free") +
    theme_linedraw() +
    plot_theme() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), labels=scales::percent_format()) +
    ylab("Relative abundance") +
    labs(fill=str_extract(plot_name, "Class|Genus")) +
    scale_fill_manual(values=make_colour_map(plot_data$name)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid=element_line(colour="black", size=0.5)) +
    xlab("Sample")
}
stacked_barplots$Class.fungal


#
# RF two-sample test
auc_plot_titles <- c(
  Class.fungal="Predicting fungal cluster from bacterial community",
  Genus.bacterial="Predicting bacterial cluster from fungal community"
)

aucs_w_pvalues <- arrow::read_feather(
  "../results/manuscript_data_finalfinal/performance_aucs_pvalues.feather"
) %>%
  right_join(
    tribble(
      ~phenotype, ~kingdom,
      "Class.fungal", "bacterial",
      "Genus.bacterial", "fungal"
    )
  ) %>%
  dplyr::rename(type=curvetypes) %>%
  group_by(phenotype, transform, k, type) %>%
  mutate(p_value=p.adjust(p_value, method="fdr")) %>%
  ungroup() %>%
  mutate(kingdom=factor(kingdom, levels=c("bacterial", "fungal")))

tick_font_size <- 10
axis_title_font_size <- 12

auroc_panels <- aucs_w_pvalues %>%
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
      scale_fill_manual(drop=FALSE, values=scales::hue_pal()(3)[1:2])
  )
auroc_panels <- setNames(auroc_panels$plot, auroc_panels$phenotype)
auroc_panels$Class.fungal

# AU PRC panels
sample_groups <- readRDS("../Data/formatted/sample_groups.rds")
prc_baseline_values <- lapply(sample_groups, function(x) as.factor(x) %>% droplevels() %>% as.integer()-1)
prc_baseline_values <- sapply(prc_baseline_values, function(x) sum(x)/length(x)) 

auprc_panels <- aucs_w_pvalues %>%
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
      scale_fill_manual(drop=FALSE, values=scales::hue_pal()(3)[1:2]) +
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
      )
  )
auprc_panels$plot[[1]]
auprc_panels <- setNames(auprc_panels$plot, auprc_panels$phenotype)

ggpubr::ggarrange(
  auroc_panels$Genus.bacterial, auprc_panels$Genus.bacterial,
  auroc_panels$Class.fungal, auprc_panels$Class.fungal,
  nrow=2, ncol=2
)

fig <- ggpubr::ggarrange(
  ggarrange(
    ggarrange(
      plotlist=lapply(biplots$plot, function(x) x + theme(title=element_blank(), legend.position="top")),
      ncol=1,
      common.legend=FALSE, 
      labels=c("a", ""),
      align="v"
    ),
    ggarrange(
      auroc_panels$Genus.bacterial, auprc_panels$Genus.bacterial,
      auroc_panels$Class.fungal, auprc_panels$Class.fungal,
      nrow=2, ncol=2,
      common.legend=TRUE, legend="right",
      labels=c("b", "", "c", "")
    ),
    ncol=2, nrow=1,
    widths=c(1, 2.5)
  ),
  ggarrange(
    stacked_barplots$Genus.bacterial,
    stacked_barplots$Class.fungal,
    # align="hv", axis="trbl",
    ncol=1,
    labels=c("d", "e")
  ),
  nrow=2, ncol=1
  # heights=c(1.5, 1, 1),
  # label.args=list(gp=grid::gpar(cex=1.2, fontface="bold"))
)

ggsave("../plots/dir_mix_figure.png", plot=fig, height=10, width=12, dpi=450)
