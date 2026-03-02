### At Time Course - Camalexin
### February 2026 AJM

library(tidyverse)

cam_df <- read.csv("data/genes_of_interest/cam_genes.csv")
cam <- cam_df %>% pull(gene)
df <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_adjusted_emmeans_log.csv")

#remove `.1` from proten/gene ID, e.g.
df <- df %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))

#filter to cam genes
df <- df %>%
	filter(protein_ID %in% cam)

df <- left_join(df, cam_df, join_by(protein_ID == gene))

df_centered <- df %>%
	group_by(protein_ID) %>%
	mutate(emmean_centered = emmean_log - mean(emmean_log, na.rm = TRUE)) %>%
	ungroup()

df_centered %>%
	ggplot(aes(x = hpi,
						 y = emmean_centered,
						 color = treatment,
						 group = treatment)) +
	geom_line(size = 1) +
	geom_point(size = 2) +
	facet_wrap(~ gene_name, scales = "free_y") +
	theme_minimal() +
	labs(
		x = "Time (hpi)",
		y = "Protein Abundance",
		color = "Treatment"
	) +
	theme(
		legend.position = c(0.63, 0.12),   # x, y
		plot.title = element_text(hjust = 0.5, face = "bold")
		#legend.key.size = unit(0.4, "cm")
	) +
	ggtitle("Camalexin Biosynthetic Proteins") +
	scale_color_brewer(palette = "Dark2", direction = -1)

ggsave("figures/timecourse/arabidopsis/metabolites/cam.png",
			 height = 4.5, width = 7.5)
