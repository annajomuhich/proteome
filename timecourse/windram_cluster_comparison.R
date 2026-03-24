### Comparing At Protein Time Course to Windram clusters
### February 2026

library(tidyverse)

### ---- Which clusters do our proteins belong to?  -------
clust <- read.csv("data/Windrametal2012/windram_clusters.csv")
#genes are unique and non-overlapping across clusters
unique(clust$AGI) %>% length()

anova <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_anova.csv")
#remove `.1` from proten/gene ID, e.g.
anova <- anova %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))
prot_int <- anova %>%
	filter(variable == "treatment" | variable == "treatment:hpi") %>%
	filter(p_adj < 0.05) %>%
	pull(protein_ID) %>%
	unique()

df <- data.frame(
	protein = prot_int
)
df <- left_join(df, clust, join_by("protein" == "AGI"))

df %>%
	mutate(Cluster = ifelse(is.na(Cluster), "No cluster", as.character(Cluster))) %>%
	count(Cluster) %>%
	mutate(Cluster = reorder(Cluster, -n)) %>%
	ggplot(aes(x = Cluster, y = n)) +
	geom_col(fill = "dark blue") +
	theme_minimal() +
	labs(x = "Windram et al Cluster ID", y = "Number of Proteins") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("figures/timecourse/arabidopsis/windram_cluster_comparison/protein_distribution.png",
# 			 height = 4, width = 8)

### ------- How does each cluster profile correspond with our protein ANOVA? ------------
protein_found <- df %>%
	filter(!is.na(Cluster)) %>%
	pull(protein)

#remove single gene clusters
clust <- clust %>%
	filter(Cluster != 35) %>%
	filter(Cluster != 7)

#get significance IDs for each protein
sig_proteins <- anova %>%
	filter(p_adj < 0.05) %>%
	group_by(protein_ID) %>%
	summarise(
		significance = sort(unique(variable)) |> 
			paste(collapse = "_"),
		.groups = "drop"
	)

#join significance data with all proteins
all_proteins <- anova %>%
	pull(protein_ID) %>%
	unique()
all_proteins <- data.frame(
	protein_ID = all_proteins
)
all_proteins <- left_join(all_proteins, sig_proteins, by = "protein_ID")
#convert NAs to "not significant"
all_proteins <- all_proteins %>%
	mutate(significance = if_else(is.na(significance), "not significant", significance))

#join this info to the cluster data
clust <- left_join(clust, all_proteins, join_by("AGI" == "protein_ID"))

#convert NAs to "not detected"
clust <- clust %>%
	mutate(significance = if_else(is.na(significance), "not detected", significance))

#reorder the significance fctors
clust <- clust %>%
	mutate(significance = factor(significance,
															 levels = c("not detected",
															 					 "not significant",
															 					 "hpi",
															 					 "treatment",
															 					 "treatment:hpi",
															 					 "hpi_treatment",
															 					 "hpi_treatment:hpi",
															 					 "treatment_treatment:hpi",
															 					 "hpi_treatment_treatment:hpi")))

clust %>%
	mutate(Cluster = as.factor(Cluster),
				 Cluster = forcats::fct_rev(Cluster)) %>%
	count(Cluster, significance) %>%
	ggplot(aes(x = n, y = Cluster, fill = significance)) +
	geom_col() +
	theme_minimal() +
	labs(
		x = "Gene count",
		y = "Windram et al Cluster",
		fill = "Protein significance"
	) +
	scale_fill_brewer(
		palette = "Spectral",
		labels = function(x) {
			x <- gsub("hpi", "time", x)
			x <- gsub("treatment", "infection", x)
			gsub("_", ", ", x)
		}
	) +
	theme(
		legend.position = c(0.75, 0.27),   # x, y
		legend.background = element_rect(fill = "white", color = "black"),
		legend.key.size = unit(0.4, "cm")
	)
ggsave("figures/timecourse/arabidopsis/windram_cluster_comparison/total_cluster_membership_sig.png",
			 width = 5, height = 5)

clust <- clust %>%
	mutate(protein_found = if_else(AGI %in% protein_found, "yes", "no"))

# plot_df <- clust %>%
# 	mutate(
# 		Cluster = as.character(Cluster),
# 		protein_found = factor(protein_found, levels = c("yes", "no"))
# 	) %>%
# 	count(Cluster, protein_found) %>%
# 	group_by(Cluster) %>%
# 	mutate(
# 		total = sum(n),
# 		prop_yes = n[protein_found == "yes"] / total
# 	) %>%
# 	ungroup() %>%
# 	mutate(Cluster = reorder(Cluster, -prop_yes))
# 
# ggplot(plot_df, aes(x = Cluster, y = n, fill = protein_found)) +
# 	geom_col() +
# 	geom_text(
# 		data = distinct(plot_df, Cluster, total, prop_yes),
# 		aes(x = Cluster, y = total, 
# 				label = scales::percent(prop_yes, accuracy = 1)),
# 		vjust = -0.3, size = 2,
# 		inherit.aes = FALSE
# 	) +
# 	theme_minimal() +
# 	labs(
# 		x = "Cluster ID",
# 		y = "Number of Genes in Cluster",
# 		fill = "Protein detected"
# 	) +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# 	scale_fill_manual(
# 		values = c(
# 			"yes" = "dark blue",   # blue
# 			"no"  = "dark grey"    # red
# 		)
# 	) 
# ggsave("figures/timecourse/arabidopsis/windram_cluster_comparison/total_cluster_membership.png",
# 			 width = 8, height = 4)

### -------- Protein in cluster expression over time -----------------------------------------

#prepare protein abundance data
emmeans <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_adjusted_emmeans_log.csv")
emmeans <- emmeans %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))
emmeans <- emmeans %>%
	filter(treatment == "infected")
emmeans$hpi <- as.factor(emmeans$hpi)

#get list of clusters to loop through
clusters <- clust$Cluster %>% unique()

#initialize list to put plots in
plots <- list()

for (cluster in clusters) {
  prot_inclust <- clust %>%
    filter(Cluster == cluster) %>%
    filter(protein_found == "yes") %>%
    pull(AGI)
  plot_data <- emmeans %>%
    filter(protein_ID %in% prot_inclust) %>%
    group_by(protein_ID) %>%
    mutate(centered_emmean_log = emmean_log - mean(emmean_log, na.rm = TRUE)) %>%
    ungroup()
  
  summary_data <- plot_data %>%
    group_by(hpi) %>%
    summarise(
      mean_emmean = mean(centered_emmean_log, na.rm = TRUE),
      sd_emmean   = sd(centered_emmean_log, na.rm = TRUE),
      .groups = "drop"
    )
  
  plots[[length(plots) + 1]] <- ggplot(plot_data,
                                       aes(x = hpi,
                                           y = centered_emmean_log,
                                           group = protein_ID)) +
    # individual protein lines
    geom_line(color = "grey80", alpha = 0.6) +
    # mean line
    geom_line(data = summary_data,
              aes(x = hpi, y = mean_emmean, group = 1),
              color = "blue",
              linewidth = 1.2) +
    # + SD
    geom_line(data = summary_data,
              aes(x = hpi, y = mean_emmean + sd_emmean, group = 1),
              color = "blue",
              linetype = "dashed") +
    # - SD
    geom_line(data = summary_data,
              aes(x = hpi, y = mean_emmean - sd_emmean, group = 1),
              color = "blue",
              linetype = "dashed") +
    theme_minimal() +
    ylab("Centered emmean (log)") +
    xlab("Time (HPI)") +
    ggtitle(paste0("Cluster ", cluster))
}
pdf("figures/timecourse/arabidopsis/windram_cluster_comparison/protein_cluster_emmeans_meanscentered.pdf",
		width = 4, height = 2)
for (p in plots) {
	print(p)
}
dev.off()
