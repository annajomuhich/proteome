### Heatmap with HCA for Botrytis proteins on Arabidopsis mutants
### March 2026 AJM

library(tidyverse)
library(pheatmap)
library(dynamicTreeCut)
library(Polychrome)

df <- read.csv("data/mutant/input/Bc_Proteome_Mut_filtered.csv") %>%
	filter(treatment == "infected")
df_means <- df %>%
	group_by(gene_ID, genotype) %>%
	summarise(abundance_mean = mean(abundance))
df_means <- df_means %>%
	mutate(genotype = as.factor(genotype))

#z-scale within protein_ID
df_means <- df_means %>%
	group_by(gene_ID) %>%
	mutate(abundance_mean_z = as.numeric(scale(abundance_mean))) %>%
	ungroup()

#convert to matrix
heat_mat <- df_means %>%
	dplyr::select(-abundance_mean) %>%
	pivot_wider(names_from = genotype, values_from = abundance_mean_z) %>%
	column_to_rownames("gene_ID") %>%
	dplyr::select(Col0, Ler, Sha, AOP2, myb28, myb29, myb2829, tgg12) %>%
	as.matrix()

protein_dist  <- dist(heat_mat)
protein_clust <- hclust(protein_dist, method = "complete")

### Cut tree into clusters -----------------------

clusters <- cutreeDynamic(
	dendro = protein_clust,
	distM = as.matrix(dist(heat_mat)),
	#minClusterSize = 5,
	deepSplit = 1,
	cutHeight = 10
)

cluster_df <- data.frame(
	protein_ID = rownames(heat_mat),
	cluster = clusters
)

# order rows the same way as the heatmap
ord <- protein_clust$order
clusters_ord <- clusters[ord]

annotation_row <- data.frame(cluster = factor(clusters))
rownames(annotation_row) <- rownames(heat_mat)

set.seed(407)
cluster_colors <- createPalette(n_distinct(annotation_row$cluster), seedcolors = c("#5F4690","#1D6996","#38A6A5"))

annotation_colors <- list(
	cluster = setNames(cluster_colors, sort(unique(annotation_row$cluster)))
)

colnames(heat_mat)

png(filename = "figures/mutant/botrytis/tc_bc_heatmap_complete.png",
		width = 5000, height = 8000, res = 1500)
pheatmap(
	heat_mat,
	cluster_rows = protein_clust,
	cluster_cols = FALSE,
	scale = "none",   # <- important
	#gaps_row = gaps,
	annotation_row = annotation_row,
	annotation_colors = annotation_colors,
	show_rownames = FALSE,
	color = colorRampPalette(
		c("blue", "white", "red")
	)(50),
	#breaks = seq(-4, 4, length.out = 51),
	legend = TRUE,
	labels_col = c("Col0", "Ler", "Sha", "AOP2", "myb28", "myb29", "myb28/29", "tgg1/2")
)
dev.off()

