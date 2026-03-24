### Making a supplementary table for the model output - Arabidopsis TC
### March 2026 AJM

emmeans <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_adjusted_emmeans_response.csv")
anova <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_anova.csv")
annot <- read_tsv("data/gene_descriptions/Ath_gene_descriptions.txt")
clust <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_clusters.csv") %>%
	dplyr::select(!color)

emmeans <- emmeans %>%
	mutate(sample_ID = paste0(treatment, "_", hpi, "HPI_mean")) %>%
	dplyr::select(emmean_response, sample_ID, protein_ID) %>%
	pivot_wider(names_from = sample_ID, values_from = emmean_response)

anova <- anova %>%
    dplyr::select(protein_ID, variable, p_adj) %>%
    pivot_wider(names_from = variable, values_from = p_adj) %>%
    dplyr::rename(pval_time = hpi,
    pval_infection = treatment,
    pval_interaction = `treatment:hpi`) %>%
	dplyr::select(!intercept)

annot <- annot %>%
	dplyr::select(`Gene Model Name`, `Primary Gene Symbol`)

df <- left_join(anova, annot, join_by(protein_ID == `Gene Model Name`))

df <- left_join(df, clust, by = "protein_ID")

df <- left_join(df, emmeans, by = "protein_ID")

df <- df %>%
	dplyr::select(protein_ID,
								`Primary Gene Symbol`,
								cluster, everything()) %>%
	dplyr::rename(gene_description = `Primary Gene Symbol`)

df %>%
	write.csv("tables/tc_at_model_table.csv", row.names = F)
