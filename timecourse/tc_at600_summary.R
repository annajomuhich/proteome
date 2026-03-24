### ANOVA summary of ~600 genes that were modeled within timepoint
### March 2026 AJM

library(tidyverse)

df <- read.csv("data/timecourse/tc_at_model_timepoint_20260302/at_tc_anova.csv")

df <- df %>%
	filter(variable == "treatment",
				 convergence_note == "relative convergence (4)") 

# Summarize protein counts for each hpi
summary_df <- df %>%
  group_by(hpi) %>%
  summarise(
    n_total = n_distinct(protein_ID),
    n_signif = n_distinct(protein_ID[p_adj < 0.05])
  ) %>%
  mutate(
    n_non_signif = n_total - n_signif
  )

# Prepare data for stacked bar plot
plot_df <- summary_df %>%
  dplyr::select(hpi, n_signif, n_non_signif) %>%
  tidyr::pivot_longer(
    cols = c(n_signif, n_non_signif),
    names_to = "category",
    values_to = "count"
  ) %>%
  mutate(
    category = factor(category, levels = c("n_non_signif", "n_signif"),
                     labels = c("Not Significant", "Significant"))
  )

ggplot(
  plot_df %>% mutate(category = forcats::fct_relevel(category, "Significant", "Not Significant")), 
  aes(x = factor(hpi), y = count, fill = category)
) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = c("Not Significant" = "grey80", "Significant" = "#155799"),
    name = ""
  ) +
  labs(
    x = "Time (hpi)",
    y = "Number of proteins"
  ) +
  theme_minimal()
ggsave("figures/timecourse/arabidopsis/tc_at600_summary.png", height = 5, width = 5)

### -------- What are the 10 significant 48H proteins? -------------------------

annot <- read.delim("data/gene_descriptions/Ath_gene_descriptions.txt")

sig10 <- df %>%
	filter(hpi == "48",
				 p_adj < 0.05)
sig10 <- left_join(sig10, annot, join_by(protein_ID == Locus.Identifier))
sig10 <- sig10 %>%
	dplyr::select(protein_ID, p_adj, Gene.Model.Type)

sig10 %>% write.csv("data/timecourse/tc_at_model_timepoint_20260302/48hpi_sig10_annots.csv", row.names = F)
