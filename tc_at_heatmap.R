### Heatmap with HCA for the 2293 infection-responsive At proteins
### February 2026 AJM

library(tidyverse)
library(pheatmap)

df <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_adjusted_emmeans_log.csv")
prot_int <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_anova.csv") %>%
	filter(variable == "treatment" | variable == "treatment:hpi") %>%
	filter(p_adj < 0.05) %>%
	pull(protein_ID) %>%
	unique()

df <- df %>%
	filter(protein_ID %in% prot_int)

#convert to factor
df <- df %>%
	mutate(treatment = as.factor(treatment),
				 hpi = as.factor(hpi))

#compute logFC
fc_df <- df %>%
	select(protein_ID, treatment, hpi, emmean_log) %>%
	pivot_wider(names_from = treatment, values_from = emmean_log) %>%
	mutate(log2FC = (infected - mock)/log(2)) %>%
	select(protein_ID, hpi, log2FC)

#sample a subset of proteins
sampled_proteins <- fc_df %>%
	distinct(protein_ID) %>%
	slice(seq(1, n(), by = 4)) %>%
	pull(protein_ID)

fc_sampled <- fc_df %>%
	filter(protein_ID %in% sampled_proteins)

#convert to matrix
heat_mat <- fc_sampled %>%
	pivot_wider(names_from = hpi, values_from = log2FC) %>%
	column_to_rownames("protein_ID") %>%
	as.matrix()

protein_dist  <- dist(heat_mat)
protein_clust <- hclust(protein_dist, method = "complete")

png(filename = "figures/timecourse/arabidopsis/tc_at_heatmap.png",
		width = 4000, height = 8000, res = 1500)
pheatmap(
	heat_mat,
	cluster_rows = protein_clust,
	cluster_cols = FALSE,
	scale = "none",   # <- important
	show_rownames = FALSE,
	color = colorRampPalette(
		c("blue", "white", "red")
	)(50),
	breaks = seq(-4, 4, length.out = 51),
	legend = TRUE,
	labels_col = c("0", "8", "16", "24", "32", "40", "48")
)
dev.off()

### --------------- pull out 48hpi upreg cluster ---------------------------

protein_order <- protein_clust$labels[protein_clust$order]
bottom_block <- tail(protein_order, 120) #this pulls the bottom ~quarter of the heatmap
#value is a little subjective but whatever

#yes, these have high values at 48HAI and look right
heat_mat_check <- heat_mat %>%
	as.data.frame() %>%
	rownames_to_column(var = "protein_ID") %>%
	filter(protein_ID %in% bottom_block)

#saving this cluster for GO analysis
high_48hpi_cluster <- heat_mat_check %>% pull(protein_ID) %>% as.data.frame()
colnames(high_48hpi_cluster) <- "protein_ID"
high_48hpi_cluster %>% write.csv("data/genes_of_interest/high_48hpi_cluster.csv", row.names = F)
