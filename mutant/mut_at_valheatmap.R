### Arabidopsis Mutant Proteome Validation Heatmap
### March 2026 AJM

library(tidyverse)

df <- read.csv("data/mutant/input/At_Proteome_Mut_filtered.csv")
gsl <- read.csv("data/genes_of_interest/gsl_forheatmap.csv")

gsl_genes <- gsl %>% pull(AGI)
df <- df %>%
	filter(gene_ID %in% gsl_genes)

df <- df %>%
    group_by(gene_ID, genotype, treatment) %>%
    summarise(abundance_mean = mean(abundance)) %>%
    mutate(genotype = factor(genotype, levels = c("Col0", "Ler", "Sha", "AOP2", "myb28", "myb29", "myb2829", "tgg12")))

df <- left_join(df, gsl, join_by(gene_ID == AGI))

# library(pheatmap)
# 
# # Filter to only infected treatment
# heat_mat <- df %>%
#   filter(treatment == "infected") %>%
#   group_by(GeneName) %>%
#   mutate(abundance_mean_z = as.numeric(scale(abundance_mean))) %>%
#   ungroup() %>%
#   mutate(genotype = as.character(genotype)) %>%
#   dplyr::select(GeneName, genotype, abundance_mean_z) %>%
#   pivot_wider(names_from = genotype, values_from = abundance_mean_z) %>%
#   column_to_rownames("GeneName") %>%
#   dplyr::select(Col0, Ler, Sha, AOP2, myb28, myb29, myb2829, tgg12) %>%
#   as.matrix()
# 
# # Convert NAs to 0
# heat_mat[is.na(heat_mat)] <- 0
# 
# # Reorder the row names of heat_mat to match the order of GeneName in gsl
# gsl_names_order <- gsl %>% pull(GeneName)
# heat_mat <- heat_mat[gsl_names_order[gsl_names_order %in% rownames(heat_mat)], , drop = FALSE]
# 
# 
# # Plot the heatmap with customized column labels
# pheatmap(
#   heat_mat,
#   cluster_rows = FALSE,
#   cluster_cols = FALSE,
#   scale = "none",
#   color = colorRampPalette(c("blue", "white", "red"))(50),
#   show_rownames = TRUE,
#   show_colnames = TRUE,
#   main = "Infected",
#   labels_col = c("Col0", "Ler", "Sha", "AOP2", "myb28", "myb29", "myb28/29", "tgg1/2")
# )

heat_df <- df %>%
	group_by(GeneName) %>%
	mutate(abundance_mean = as.numeric(scale(abundance_mean))) %>%
	ungroup() %>%
	mutate(genotype = as.character(genotype)) %>%
	dplyr::select(treatment, GeneName, genotype, abundance_mean)

# enforce ordering
gsl_names_order <- gsl %>% pull(GeneName)

heat_df <- heat_df %>%
	filter(GeneName %in% gsl_names_order) %>%
	mutate(GeneName = factor(GeneName, levels = gsl_names_order),
				 genotype = factor(
				 	genotype, 
				 	levels =c("Col0", "Ler", "Sha", "AOP2", "myb28", "myb29", "myb2829", "tgg12")),
				 treatment = factor(
				 	treatment,
				 	levels = c("mock", "infected")))

ggplot(heat_df, aes(x = genotype, y = GeneName, fill = abundance_mean)) +
	geom_tile() +
	facet_wrap(~ treatment) +
	scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
	scale_x_discrete(labels = c("Col0", "Ler", "Sha", "AOP2", "myb28", "myb29", "myb28/29", "tgg1/2")) +
	scale_y_discrete(limits = rev) +
	labs(x = "", y = "", fill = "Z-score") +
	theme_minimal() +
	theme(
		axis.text.y = element_text(size = 8),
		axis.text.x = element_text(angle = 45, hjust = 1),
		strip.text = element_text(face = "bold")
	)
ggsave("figures/mutant/arabidopsis/mut_at_valheatmap.png", height = 5.16, width = 4.3,
			dpi = 1500)

#presence/absence option
#do before scaling the data
# 
# heat_df_pa <- heat_df %>%
# 	mutate(presence = ifelse(abundance_mean != 0 & !is.na(abundance_mean), 1, 0))
# ggplot(heat_df_pa, aes(x = genotype, y = GeneName, fill = factor(presence))) +
# 	geom_tile(color = "grey90") +
# 	facet_wrap(~ treatment) +
# 	scale_fill_manual(
# 		values = c("0" = "white", "1" = "black"),
# 		labels = c("Absent", "Present"),
# 		name = ""
# 	) +
# 	scale_x_discrete(labels = c("Col0", "Ler", "Sha", "AOP2", "myb28", "myb29", "myb28/29", "tgg1/2")) +
# 	scale_y_discrete(limits = rev) +
# 	labs(x = "", y = "") +
# 	theme_minimal() +
# 	theme(
# 		axis.text.y = element_text(size = 8),
# 		axis.text.x = element_text(angle = 45, hjust = 1),
# 		strip.text = element_text(face = "bold")
# 	)
