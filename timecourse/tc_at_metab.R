### Plotting aliphatic vs indole proteins colored by cluster
### March 2026 AJM

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
