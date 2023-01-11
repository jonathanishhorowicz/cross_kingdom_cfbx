library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(tidybayes)
library(stringr)
library(dendextend)
library(ComplexHeatmap)
library(corrr)
library(ggplot2)
library(SpiecEasi)

options(stringsAsFactors=FALSE)

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

taxonomy <- lapply(
  phylos,
  function(x) tax_table(x)@.Data %>% as_tibble()
) %>%
  bind_rows() 

# agglomerate to Genus
cat(sprintf("Agglomerating to %s\n", "Genus"))
phylos_glom <- lapply(phylos, function(x) speedyseq::tax_glom(x, "Genus"))
sapply(phylos_glom, ntaxa)

# get OTU table and remove zero-variance OTUs
X <- lapply(
  phylos_glom,
  function(x) x %>% otu_table() %>% as.data.frame() %>% t()
) 
X <- lapply(X, function(x) x[ rownames(y), ])
X <- lapply(X, function(x) x[,matrixStats::colVars(x)!=0])

# prevalence filtering - remove OTUs present in fewer than 10% of samples
taxa_prevalences <- lapply(X, function(x) apply(x, 2, function(x) sum(x>0)))
prevalence_cutoff <- 10
included_taxa <- lapply(X, function(x) colnames(x)[apply(x, 2, function(xx) sum(xx>0))>prevalence_cutoff]) 
lapply(included_taxa, length)

all_included_taxa <- lapply(included_taxa, as_tibble) %>%
  bind_rows(.id="kingdom") %>%
  rename(taxa=value) %>%
  left_join(taxonomy, by=c(taxa="OTUID"))

# store the kingdom and plotted name of each taxa
taxa_renamer <- setNames(str_remove(all_included_taxa$Genus, "f__|g__|o__|c__"), all_included_taxa$taxa)
kingdom_getter <- taxonomy %>%
  select(Genus, Domain) %>%
  distinct()
kingdom_getter <- setNames(kingdom_getter$Domain, str_remove(kingdom_getter$Genus, "f__|g__|o__|c__"))

###############################################################
###############################################################

###############################################################
# run sparCC on the agglomerated, CLR-transformed counts
###############################################################

X_sparcc <- do.call(cbind, X)[,all_included_taxa$taxa]
dim(X_sparcc)

# sparCC correlations
sparcc_res <- sparcc(X_sparcc)

###############################################################
###############################################################

###############################################################
# make heatmaps
###############################################################

# only plot high prevalence, high abundance taxa in main text heatmap
main_text_plot_cutoff <- 25
all_taxa_prevalence <- c(taxa_prevalences$fungal, taxa_prevalences$bacterial)
plotted_taxa_maintext <-  names(all_taxa_prevalence)[ all_taxa_prevalence > main_text_plot_cutoff ]


# heatmap image
im <- sparcc_res$Cor
rownames(im) <- colnames(X_sparcc)
colnames(im) <- colnames(X_sparcc)
im <- im[ plotted_taxa_maintext,plotted_taxa_maintext ]

# rename using genus names
rownames(im) <- taxa_renamer[ rownames(im) ]
colnames(im) <- taxa_renamer[ colnames(im) ]

# create ComplexHeatmap, starting with annotations
annotation_data <- as_tibble(plotted_taxa_maintext) %>%
  rename(taxa=value) %>%
  left_join(taxonomy, by=c(taxa="OTUID")) %>%
  mutate(taxa=taxa_renamer[ taxa ]) %>%
  arrange(match(taxa, rownames(im))) %>%
  mutate(kingdom=kingdom_getter[ taxa ] %>%
           str_remove("^k__")) %>%
  select(kingdom, taxa)

# kingdom annotation colour palette
kingdom_cmap <- c(
  Fungi="#00BA38",
  Bacteria="#F8766D"
)

# sparCC correlations heatmap (panel a)
hm <- Heatmap(
  im,
  name="SparCC",
  col=circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  heatmap_height=unit(8, "inches"),
  heatmap_width=unit(8, "inches"),
  top_annotation=columnAnnotation(kingdom=annotation_data$kingdom,
                                  col=list(kingdom=kingdom_cmap),
                                  show_legend=FALSE),
  left_annotation=rowAnnotation(kingdom=annotation_data$kingdom,
                                col=list(kingdom=kingdom_cmap)),
  row_names_gp=gpar(fontsize=5),
  column_names_gp=gpar(fontsize=5),
  row_dend_width=unit(30, "mm"),
  column_dend_height=unit(30, "mm")
)

pdf("../plots/sparcc_heatmap_supplemental.pdf", width=18, height=15)
draw(hm, merge_legend=TRUE)
dev.off()
knitr::plot_crop("../plots/sparcc_heatmap_supplemental.pdf")

hm_grob <- grid::grid.grabExpr(draw(hm, merge_legend=TRUE, padding=unit(c(-10, 0, 10, 0), "mm")))

# sparCC p-values plot (panel b)
im2 <- sparcc_res$Cor
rownames(im2) <- taxa_renamer[ colnames(X_sparcc) ]
colnames(im2) <- taxa_renamer[ colnames(X_sparcc) ]

plt_b <- abs(im2) %>%
  as_cordf() %>%
  shave() %>%
  stretch() %>%
  mutate(
    facet_title=sprintf("%s-%s", kingdom_getter[ x ], kingdom_getter[ y ]) %>%
      str_remove_all("k__") %>%
      str_replace("-", "\n") %>%
      factor(levels=c("Bacteria\nBacteria", "Fungi\nFungi", "Fungi\nBacteria"))
  ) %>%
  filter(!is.na(r)) %>%
  ggplot(aes(x=facet_title, y=r)) +
  geom_boxplot() +
  theme_linedraw() +
  theme(
    axis.title.x=element_blank(),
    axis.text=element_text(size=10),
    axis.title.y=element_text(size=14),
    panel.grid.minor.y=element_line(size=0.005, colour="grey50"),
    panel.grid.major.y=element_line(size=0.005, colour="grey50"),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank()
  ) +
  ylab("|SparCC correlation|")