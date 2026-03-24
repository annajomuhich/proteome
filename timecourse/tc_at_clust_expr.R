### Arabidopsis time course: Plotting cluster expression over time
### March 2026 AJM

library(dplyr)
library(ggplot2)

clust <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_clusters.csv")
expr <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_adjusted_emmeans_response.csv")

# Convert cluster to character to maintain consistency
clust$cluster <- as.character(clust$cluster)

expr_scaled <- expr %>%
	group_by(protein_ID) %>%
	mutate(emmean_scaled = scale(emmean_response),
				 hpi = as.factor(hpi))

for (cl in unique(clust$cluster)) {
	# Pull protein_IDs for this cluster
	clust_proteins <- clust %>%
		filter(cluster == cl) %>%
		pull(protein_ID)
	
	# Filter expr for these protein_IDs
	expr_clust <- expr_scaled %>%
		filter(protein_ID %in% clust_proteins)
	
	# Get color for this cluster (assumes color field exists in `clust`)
	clust_color <- clust %>% filter(cluster == cl) %>% dplyr::slice(1) %>% pull(color)
	if (length(clust_color) == 0 || is.na(clust_color)) {
		clust_color <- "#FF0000" # fallback color
	}
	
	# Calculate mean and standard deviation for each hpi and treatment
	expr_mean <- expr_clust %>%
		group_by(hpi, treatment) %>%
		summarize(
			mean_emmean = mean(emmean_scaled, na.rm = TRUE),
			sd_emmean = sd(emmean_scaled, na.rm = TRUE),
			.groups = "drop"
		)
	
	# Plot just the mean, with a grey ribbon for +/- 1 sd
	p <- ggplot(expr_mean, aes(x = hpi, y = mean_emmean, group = treatment, color = treatment, linetype = treatment)) +
		geom_ribbon(aes(ymin = mean_emmean - sd_emmean, ymax = mean_emmean + sd_emmean, fill = treatment), 
								alpha = 0.2, color = NA, show.legend = FALSE) +
		geom_line(size = 1.3) +
		scale_linetype_manual(
			values = c("infected" = "solid", "mock" = "dashed"),
			breaks = c("infected", "mock")
		) +
		scale_color_manual(
			values = c("infected" = clust_color, "mock" = clust_color),
			guide = "none"
		) +
		scale_fill_manual(
			values = c("infected" = "grey80", "mock" = "grey80"),
			guide = "none"
		) +
		theme_minimal() +
		labs(
			title = paste("Cluster", cl),
			x = NULL,
			y = NULL
		) +
		theme(
			legend.position = "none"
		)
	# Save (optional, remove if not desired)
	ggsave(
		sprintf("figures/timecourse/arabidopsis/cluster_expr/cluster_%s_expression.png", cl),
		p, width = 3.0, height = 2.25
	)
}

# ### Plot all genes in each cluster (messy)

# for (cl in unique(clust$cluster)) {
#   # Pull protein_IDs for this cluster
#   clust_proteins <- clust %>%
#     filter(cluster == cl) %>%
#     pull(protein_ID)
#   
#   # Filter expr for these protein_IDs
#   expr_clust <- expr %>%
#     filter(protein_ID %in% clust_proteins)
#   
#   # Get color for this cluster (assumes color field exists in `clust`)
#   clust_color <- clust %>% filter(cluster == cl) %>% dplyr::slice(1) %>% pull(color)
#   if (length(clust_color) == 0 || is.na(clust_color)) {
#     clust_color <- "#FF0000" # fallback color
#   }
#   
#   # Calculate mean response for each hpi and treatment
#   expr_mean <- expr_clust %>%
#     group_by(hpi, treatment) %>%
#     summarize(mean_emmean = mean(emmean_log, na.rm = TRUE), .groups = "drop")
#   
#   # Plot
#   p <- ggplot(expr_clust, aes(x = hpi, y = emmean_log, group = interaction(protein_ID, treatment))) +
#     geom_line(aes(linetype = treatment), color = "grey80", size = 0.5, alpha = 0.75) +
#     geom_line(
#       data = expr_mean,
#       aes(x = hpi, y = mean_emmean, group = treatment, color = treatment, linetype = treatment),
#       size = 1.3
#     ) +
#     scale_linetype_manual(
#       values = c("infected" = "solid", "mock" = "dashed"),
#       breaks = c("infected", "mock")
#     ) +
#     scale_color_manual(
#       values = c("infected" = clust_color, "mock" = clust_color),
#       guide = "none"
#     ) +
#     theme_minimal() +
#     labs(
#       title = paste("Cluster", cl, "Expression"),
#       x = "Hours post infection (hpi)",
#       y = "Expression (emmean_response)"
#     ) +
#     theme(
#       legend.position = "none"
#     )
# 
#   # Save (optional, remove if not desired)
#   ggsave(
#     sprintf("figures/timecourse/arabidopsis/cluster_expr/cluster_%s_expression.png", cl),
#     p, width = 6, height = 4
#   )
# }
