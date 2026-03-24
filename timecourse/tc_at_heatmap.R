### Heatmap with HCA for the infection-responsive At proteins
### February 2026 AJM

#remotes::install_version("Polychrome", version = "1.5.1")
# BiocManager::install(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"))
# install.packages("WGCNA")
library(tidyverse)
library(pheatmap)
library(dynamicTreeCut)
library(Polychrome)

df <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_adjusted_emmeans_log.csv")

#Get list of significant proteins
prot_int <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_anova.csv") %>%
	filter(variable == "treatment" | variable == "treatment:hpi") %>%
	filter(p_adj < 0.05) %>%
	pull(protein_ID) %>%
	unique()

# #Collect 1 timepoint only proteins. There are 957 of them
# prot_1tp <- read.csv('data/timecourse/input/AtBc_Proteome_TimeCourse_filtered.csv') %>%
# 	filter(abundance >= 0.8) %>%
# 	group_by(protein_ID, treatment) %>%
# 	summarise(n_hpi = n_distinct(hpi), .groups = "drop") %>%
# 	filter(n_hpi == 7) %>%
# 	pull(protein_ID) %>%
# 	unique()

#filter to significant proteins
df <- df %>%
	filter(protein_ID %in% prot_int)
# # #filter to remove 1tp only proteins
# df <- df %>%
# 	filter(!(protein_ID %in% prot_1tp))

#convert to factor
df <- df %>%
	mutate(treatment = as.factor(treatment),
				 hpi = as.factor(hpi))

#compute logFC
fc_df <- df %>%
	dplyr::select(protein_ID, treatment, hpi, emmean_log) %>%
	pivot_wider(names_from = treatment, values_from = emmean_log) %>%
	mutate(log2FC = (infected - mock)/log(2)) %>% #this one if values are already logged
	#mutate(log2FC = log2(infected/mock)) #%>%		#this one if values are response-scale
	dplyr::select(protein_ID, hpi, log2FC)

# Convert outlier log2FC (from presence/absence expression)
# Log2FC values will max out at 5 and min at -5 for visualization
fc_df <- fc_df %>%
	mutate(log2FC = ifelse(log2FC > 5, 5, ifelse(log2FC < -5, -5, log2FC)))

#sample a subset of proteins
sampled_proteins <- fc_df %>%
	distinct(protein_ID) %>%
	dplyr::slice(seq(1, n(), by = 1)) %>%
	pull(protein_ID)

fc_sampled <- fc_df %>%
	filter(protein_ID %in% sampled_proteins)

#convert to matrix
heat_mat <- fc_sampled %>%
	pivot_wider(names_from = hpi, values_from = log2FC) %>%
	column_to_rownames("protein_ID") %>%
	as.matrix()

protein_dist  <- dist(heat_mat)

#remove outliers.
avg_dist <- rowMeans(as.matrix(protein_dist))
outliers <- names(avg_dist[avg_dist > quantile(avg_dist, 0.99)])
outliers <- c(outliers, "AT4G23440.1", "AT2G22970.1", "AT1G14870.1")
# #make sure camalexin genes are not removed
# outliers <- outliers[outliers != "AT3G26830.1"]
# outliers <- outliers[outliers != "AT2G30750.1"]
# outliers <- outliers[outliers != "AT2G30770.1"]

heat_mat <- heat_mat %>% 
	as.data.frame() %>% 
	rownames_to_column("protein_ID") %>%
	filter(!(protein_ID %in% outliers)) %>%
	column_to_rownames("protein_ID") %>%
	as.matrix()

protein_dist  <- dist(heat_mat)
protein_clust <- hclust(protein_dist, method = "complete")

### Cut tree into clusters -----------------------

clusters <- cutreeDynamic(
	dendro = protein_clust,
	distM = as.matrix(dist(heat_mat)),
	minClusterSize = 30,
	deepSplit = 2
	#cutHeight = 10
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

png(filename = "figures/timecourse/arabidopsis/heatmap/tc_at_heatmap_complete_clust2_30_20260317.png",
		width = 4500, height = 8000, res = 1500)
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
	labels_col = c("0", "8", "16", "24", "32", "40", "48")
)
dev.off()

### Write out cluster assignments with their color hex codes
color_df <- as.data.frame(cluster_colors) %>%
	rownames_to_column("cluster") %>%
	mutate(cluster = gsub("NC", "", cluster)) %>%
	dplyr::rename(color = cluster_colors)
cluster_df <- cluster_df %>%
	mutate(cluster = as.character(cluster)) %>%
	left_join( color_df, by = "cluster")
cluster_df %>% write.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_clusters.csv", row.names = F)

# ### --------------- pull out 48hpi upreg cluster ---------------------------
# 
# protein_order <- protein_clust$labels[protein_clust$order]
# bottom_block <- tail(protein_order, 120) #this pulls the bottom ~quarter of the heatmap
# #value is a little subjective but whatever
# 
# #yes, these have high values at 48HAI and look right
# heat_mat_check <- heat_mat %>%
# 	as.data.frame() %>%
# 	rownames_to_column(var = "protein_ID") %>%
# 	filter(protein_ID %in% bottom_block)
# 
# #saving this cluster for GO analysis
# high_48hpi_cluster <- heat_mat_check %>% pull(protein_ID) %>% as.data.frame()
# colnames(high_48hpi_cluster) <- "protein_ID"
# high_48hpi_cluster %>% write.csv("data/genes_of_interest/high_48hpi_cluster.csv", row.names = F)
