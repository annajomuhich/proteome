### Heatmap with HCA for the infection-responsive At proteins
### February 2026 AJM

library(tidyverse)
library(pheatmap)

df <- read.csv("data/timecourse/botrytis/tc_bc_model_timepoint_20260309/bc_tc_adjusted_emmeans_response.csv")

df <- df %>%
	filter(treatment == "infected") %>%
	filter(hpi != 0)

#convert to factor
df <- df %>%
	mutate(hpi = as.factor(hpi))

#z-scale within protein ID
df <- df %>%
	group_by(protein_ID) %>%
	mutate(emmean_response_z = as.numeric(scale(emmean_response))) %>%
	ungroup()

#convert to matrix
heat_mat <- df %>%
	dplyr::select(protein_ID, hpi, emmean_response_z) %>%
	pivot_wider(names_from = hpi, values_from = emmean_response_z) %>%
	column_to_rownames("protein_ID") %>%
	na.omit() %>% #remove rows with NA
	as.matrix()

protein_dist  <- dist(heat_mat)
protein_clust <- hclust(protein_dist, method = "average")

png(filename = "figures/timecourse/botrytis/tc_bc_heatmap_average_resp.png",
		width = 3000, height = 8000, res = 1400)
pheatmap(
	heat_mat,
	cluster_rows = protein_clust,
	cluster_cols = FALSE,
	scale = "none",   # <- important
	show_rownames = FALSE,
	color = colorRampPalette(
		c("blue", "white", "red")
	)(50),
	#breaks = seq(-4, 4, length.out = 51),
	legend = TRUE,
	labels_col = c("8", "16", "24", "32", "40", "48")
)
dev.off()
