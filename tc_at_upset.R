###### Arabidopsis time course - quantitative genes - Upset plot
###### AJM February 2026

library(tidyverse)
library(ComplexUpset)

##### ------------- Set up data -------------------------
df <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260309/at_tc_anova.csv") %>%
	dplyr::select(protein_ID, variable, p_adj) %>%
	filter(variable != "intercept")

any(is.na(df))
colSums(is.na(df) > 0) #see by column


df_wide <- pivot_wider(df, names_from = "variable", values_from = "p_adj")

# convert p-values to significance (TRUE/FALSE)
sig_df <- df_wide %>%
	mutate(
		`Time` = hpi < 0.05,
		`Infection` = treatment < 0.05,
		`Time X Infection` = `treatment:hpi` < 0.05
	) %>%
	dplyr::select(protein_ID, `Time`, `Infection`, `Time X Infection`)

sig_df <- sig_df %>% dplyr::rename("Time X\nInfection" = `Time X Infection`)

### Plot ============================================================
upset(
	sig_df,
	intersect = c("Time", "Infection", "Time X\nInfection"),
	annotations = list(),
	base_annotations = list(
		'Significant Arabidopsis Proteins' = intersection_size(text = list(size = 3))
	),
	width_ratio = c(0, 1),        # removes left panel entirely
	set_sizes = FALSE             # <- this removes the "Set Size" title
) +
	theme_minimal() +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		axis.text.x  = element_blank(),
		axis.ticks.x = element_blank()
	)
ggsave("figures/timecourse/arabidopsis/tc_at_upset.png", width = 4, height = 4)
