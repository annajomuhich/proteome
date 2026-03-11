###### Botrytis time course - timepoint model - time of first accumulation
###### AJM March 2026

library(tidyverse)

anova <- read.csv("data/timecourse/botrytis/tc_bc_model_timepoint_20260309/bc_tc_anova.csv") %>%
	filter(variable == "treatment")

emmeans <- read.csv("data/timecourse/botrytis/tc_bc_model_timepoint_20260309/bc_tc_adjusted_emmeans_log.csv")

first_accumulation <- anova %>%
	filter(p_adj < 0.05) %>%
	inner_join(
		emmeans %>%
			select(protein_ID, hpi, treatment, emmean_log) %>%
			pivot_wider(names_from = treatment, values_from = emmean_log) %>%
			filter(infected > mock) %>%
			select(protein_ID, hpi),
		by = c("protein_ID", "hpi")
	) %>%
	arrange(protein_ID, hpi) %>%
	group_by(protein_ID) %>%
	mutate(
		consecutive = hpi - lag(hpi) == min(diff(sort(unique(anova$hpi)))) #this detects back to back significant hpi
	) %>%
	filter(consecutive) %>%
	summarise(first_hpi = min(hpi), .groups = "drop")

first_accumulation %>%
	count(first_hpi) %>%
	mutate(first_hpi = as.factor(first_hpi)) %>%
	ggplot(aes(x = first_hpi, y = n)) +
	geom_col(fill = "black") +
	labs(
		x = "First Accumulation Time (hpi)",
		y = "Number of Botrytis proteins"
	) +
	theme_minimal()
ggsave("figures/timecourse/botrytis/tc_bc_protaccum.png", height = 4, width = 3)
