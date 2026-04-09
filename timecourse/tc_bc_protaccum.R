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
		x = "First Detection Time (HPI)",
		y = "Number of Botrytis proteins"
	) +
	theme_minimal()
ggsave("figures/timecourse/botrytis/tc_bc_protaccum.png", height = 3, width = 3)

#Save the first accumulation lists for GO analysis
xp <- read.csv("data/gene_descriptions/Bcin_XP_key.csv")
first_df <- left_join(first_accumulation, xp, join_by("protein_ID" == "XP")) %>%
	dplyr::select(gene_ID, first_hpi)
any(is.na(first_df))
first_df %>%
	write.csv("data/timecourse/botrytis/tc_bc_model_timepoint_20260309/bc_first_hpi.csv", row.names = F)

#make a version with annotations
annot <- read.csv("data/gene_descriptions/Bcin_Annotations_Full_transcript.csv") %>%
	dplyr::select(X.Gene.ID., X.Gene.Name.or.Symbol., X.PFam.Description.)
colnames(annot) <- c("gene_ID", "gene_name", "gene_desc")
first_df <- left_join(first_df, annot, by = "gene_ID")
first_df <- arrange(first_df, first_hpi)
first_df <- unique(first_df)
first_df %>%
	write.csv("data/timecourse/botrytis/tc_bc_model_timepoint_20260309/bc_first_hpi_annot.csv", row.names = F)
