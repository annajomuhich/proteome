### Plotting aliphatic vs indole proteins colored by cluster
### March 2026 AJM

library(tidyverse)

emmeans <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_adjusted_emmeans_log.csv")
cluster <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_clusters.csv")
genes <- read.csv("data/genes_of_interest/gsl_cam_genes.csv")

aliph <- genes %>%
	filter(Pathways == "Aliphatic Glucosinolate") %>%
	pull(AGI)
indole <- genes %>%
	filter(Pathways == "Indolic Glucosinolate") %>%
	pull(AGI)
cam <- genes %>%
	filter(Pathways == "Camalexin") %>%
	pull(AGI)
precursor <- genes %>%
	filter(Pathways == "Precursor") %>%
	pull(AGI)

df <- left_join(emmeans, cluster, by = "protein_ID") %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))

cluster_colors <- distinct(df, cluster, color) %>%
	mutate(cluster = ifelse(is.na(cluster), "none", cluster),
				 color = ifelse(cluster == "none", "lightgrey", color))

#aliphatic
df %>%
	# Add pathway information to df by joining with the genes data.frame
	left_join(genes %>% dplyr::select(protein_ID = AGI, Pathways), by = "protein_ID") %>%
	# Keep only proteins that are in any desired pathway
	filter(!is.na(Pathways)) %>%
	mutate(
		cluster = ifelse(is.na(cluster), "none", cluster),
		hpi = as.factor(hpi),
		cluster = as.factor(cluster)
	) %>%
	dplyr::select(treatment, hpi, emmean_log, protein_ID, cluster, color, Pathways) %>%
	pivot_wider(names_from = treatment, values_from = emmean_log) %>%
	mutate(log2FC = (infected - mock)/log(2)) %>%
	mutate(log2FC = ifelse(log2FC > 5, 5, ifelse(log2FC < -5, -5, log2FC))) %>%
	ggplot(aes(x = hpi, y = log2FC, group = protein_ID, color = cluster)) +
	geom_line() +
	scale_color_manual(
		values = setNames(
			cluster_colors$color,     # values
			cluster_colors$cluster    # names
		)
	) +
	facet_wrap(~Pathways, scales = "free_y") +
	theme_minimal()

#indolic
df %>%
	filter(protein_ID %in% indole) %>%
	mutate(
		cluster = ifelse(is.na(cluster), "none", cluster),
		hpi = as.factor(hpi),
		cluster = as.factor(cluster)
	) %>%
	dplyr::select(treatment, hpi, emmean_log, protein_ID, cluster, color) %>%
	pivot_wider(names_from = treatment, values_from = emmean_log) %>%
	mutate(log2FC = (infected - mock)/log(2)) %>%
	mutate(log2FC = ifelse(log2FC > 5, 5, ifelse(log2FC < -5, -5, log2FC))) %>%
	ggplot(aes(x = hpi, y = log2FC, group = protein_ID, color = cluster)) +
	geom_line() +
	scale_color_manual(
		values = setNames(
			cluster_colors$color,     # values
			cluster_colors$cluster    # names
		)
	) +
	theme_minimal()

#camalexin
df %>%
	filter(protein_ID %in% cam) %>%
	mutate(
		cluster = ifelse(is.na(cluster), "none", cluster),
		hpi = as.factor(hpi),
		cluster = as.factor(cluster)
	) %>%
	dplyr::select(treatment, hpi, emmean_log, protein_ID, cluster, color) %>%
	pivot_wider(names_from = treatment, values_from = emmean_log) %>%
	mutate(log2FC = (infected - mock)/log(2)) %>%
	mutate(log2FC = ifelse(log2FC > 5, 5, ifelse(log2FC < -5, -5, log2FC))) %>%
	ggplot(aes(x = hpi, y = log2FC, group = protein_ID, color = cluster)) +
	geom_line() +
	scale_color_manual(
		values = setNames(
			cluster_colors$color,     # values
			cluster_colors$cluster    # names
		)
	) +
	theme_minimal()

#precursor
df %>%
	filter(protein_ID %in% precursor) %>%
	mutate(
		cluster = ifelse(is.na(cluster), "none", cluster),
		hpi = as.factor(hpi),
		cluster = as.factor(cluster)
	) %>%
	dplyr::select(treatment, hpi, emmean_log, protein_ID, cluster, color) %>%
	pivot_wider(names_from = treatment, values_from = emmean_log) %>%
	mutate(log2FC = (infected - mock)/log(2)) %>%
	mutate(log2FC = ifelse(log2FC > 5, 5, ifelse(log2FC < -5, -5, log2FC))) %>%
	ggplot(aes(x = hpi, y = log2FC, group = protein_ID, color = cluster)) +
	geom_line() +
	scale_color_manual(
		values = setNames(
			cluster_colors$color,     # values
			cluster_colors$cluster    # names
		)
	) +
	theme_minimal()
# df_centered <- df %>%
# 	group_by(protein_ID) %>%
# 	mutate(emmean_centered = emmean_response - mean(emmean_response, na.rm = TRUE)) %>%
# 	ungroup()
# 
# df_centered %>%
# 	filter(protein_ID %in% aliph) %>%
# 	mutate(hpi = as.factor(hpi),
# 				 cluster = as.factor(cluster)) %>%
# 	ggplot(aes(x = hpi, y = emmean_centered, group = protein_ID, color = cluster)) +
# 	geom_line(aes(color = cluster, group = protein_ID)) +
# 	scale_color_manual(values = setNames(unique(df$color), unique(df$cluster)))

### Plots separated by gene - 4/9/2026 ----------------------------

### Now just need to get the p values so we can do dotted or solid lines

#prepare p values
anova <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_anova.csv") %>%
	filter(variable == "treatment" | variable == "treatment:hpi") %>%
	dplyr::select(protein_ID, variable, p_adj) %>%
	pivot_wider(names_from = variable, values_from = p_adj) %>%
	mutate(significant = ifelse(`treatment` < 0.05 | `treatment:hpi` < 0.05, "yes", "no")) %>%
	dplyr::select(protein_ID, significant) %>%
	mutate(protein_ID = sub("\\.\\d+$", "", protein_ID))

#if there are multiple p values for the same protein_ID, take the one with significant == yes
anova <- anova %>%
  group_by(protein_ID) %>%
  filter(
    # For protein_IDs that occur more than once, keep only those with significant == "yes"
    if(n() > 1) significant == "yes" else TRUE
  ) %>%
  ungroup()

# Check for duplicates in anova$protein_ID
dup_protein_ids <- anova$protein_ID[duplicated(anova$protein_ID)]
if(length(dup_protein_ids) > 0) {
  warning("Duplicate protein_IDs found in anova: ", paste(unique(dup_protein_ids), collapse = ", "))
} else {
  message("No duplicate protein_IDs in anova$protein_ID")
}


all_genes <- genes %>%
	pull(AGI) %>%
	unique()

df_plot <- df %>%
	filter(protein_ID %in% all_genes) %>%
	left_join(anova, by = "protein_ID") %>%
	left_join(genes, join_by(protein_ID == AGI)) %>%
	mutate(
		cluster = ifelse(is.na(cluster), "none", cluster),
		hpi = as.factor(hpi), 
		cluster = as.factor(cluster)
	) %>%
	dplyr::select(GeneName, Pathways, protein_ID, treatment, hpi, emmean_log,cluster, color, significant) %>%
	pivot_wider(names_from = treatment, values_from = emmean_log) %>%
	mutate(log2FC = (infected - mock)/log(2)) %>%
	mutate(log2FC = ifelse(log2FC > 5, 5, ifelse(log2FC < -5, -5, log2FC)))

# Split gene names into two halves for 2 plots (6x5 grid each)
gene_names <- unique(df_plot$GeneName)
n_per_plot <- 31

# First set of genes
genes1 <- gene_names[1:n_per_plot]
# Second set of genes
genes2 <- gene_names[(n_per_plot + 1):length(gene_names)]

# Plot 1: first 30 genes
p1 <- df_plot %>%
	filter(GeneName %in% genes1) %>% 
	ggplot(aes(
		x = hpi,
		y = log2FC,
		group = protein_ID, 
		color = cluster,
		linetype = significant
	)) +
	geom_line() +
	scale_color_manual(
		values = setNames(
			cluster_colors$color,     # values
			cluster_colors$cluster    # names
		)
	) +
	scale_linetype_manual(
		values = c("yes" = "solid", "no" = "dashed"),
		labels = c("yes" = "Significant", "no" = "Not Significant"),
		name = "Significance"
	) +
	facet_wrap(~ GeneName, scales = "free_y", ncol = 5, nrow = 7) +
	theme_classic() +
	theme(strip.background = element_blank()) +
	labs(x = "Time (HPI)", y = "Log2FC (Mock vs Infected)")

# Plot 2: next 30 genes
p2 <- df_plot %>%
	filter(GeneName %in% genes2) %>% 
	ggplot(aes(
		x = hpi,
		y = log2FC,
		group = protein_ID, 
		color = cluster,
		linetype = significant
	)) +
	geom_line() +
	scale_color_manual(
		values = setNames(
			cluster_colors$color,
			cluster_colors$cluster
		)
	) +
	scale_linetype_manual(
		values = c("yes" = "solid", "no" = "dashed"),
		labels = c("yes" = "Significant", "no" = "Not Significant"),
		name = "Significance"
	) +
	facet_wrap(~ GeneName, scales = "free_y", ncol = 5, nrow = 7) +
	theme_classic() +
	theme(strip.background = element_blank()) +
	labs(x = "Time (HPI)", y = "Log2FC (Mock vs Infected)")

# Show plots (optionally, or ggsave them if you want)
p1
ggsave("figures/timecourse/arabidopsis/metabolites/allgslgenes1.png",
			 height = 8, width = 9, dpi = 1000)
p2
ggsave("figures/timecourse/arabidopsis/metabolites/allgslgenes2.png",
			 height = 8, width = 9, dpi = 1000)

### 2026-05-04 trying it all in one big plot

df_plot %>%
	#filter(GeneName %in% genes2) %>% 
	ggplot(aes(
		x = hpi,
		y = log2FC,
		group = protein_ID, 
		color = cluster,
		linetype = significant
	)) +
	geom_line() +
	scale_color_manual(
		values = setNames(
			cluster_colors$color,
			cluster_colors$cluster
		)
	) +
	scale_linetype_manual(
		values = c("yes" = "solid", "no" = "dashed"),
		labels = c("yes" = "Significant", "no" = "Not Significant"),
		name = "Significance"
	) +
	facet_wrap(~ GeneName, scales = "free_y", ncol = 6, nrow = 11) +
	theme_classic() +
	theme(strip.background = element_blank(),
				axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
	labs(x = "Time (HPI)", y = "Log2FC (Mock vs Infected)") +
	guides(linetype = "none")
ggsave("figures/timecourse/arabidopsis/metabolites/allgslgenes.png",
			 height = 11, width = 9, dpi = 1000)
