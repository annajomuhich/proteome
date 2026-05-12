### Plotting proteins abundance of different genes of interest in Botrytis over time
### May 2026 AJM

library(tidyverse)

df <- read.csv("data/timecourse/input/AtBc_Proteome_TimeCourse_filtered.csv")

prots <- read.csv("data/genes_of_interest/bcin_vir_genes.csv")

df_filt <- df %>%
	filter(treatment == "infected",
				 organism == "botrytis",
				 protein_ID %in% prots$protein_ID) %>%
	left_join(prots, by = "protein_ID") %>%
	filter(genotype == "Col.0") %>%
	mutate(hpi = as.factor(hpi)) 

bcboa_levels <- paste0("Bcboa", 1:17)
bcboa_levels <- bcboa_levels[!(bcboa_levels %in% c("Bcboa9", "Bcboa17"))]

# Plot for Botcinic acid
df_filt %>%
  filter(category == "Botcinic acid") %>%
  filter(gene_name %in% bcboa_levels) %>%
  mutate(gene_name = factor(gene_name, levels = bcboa_levels)) %>%
  ggplot(aes(x = hpi, y = abundance, color = gene_name)) +
  geom_boxplot() +
  scale_color_brewer(palette = "Paired") +
  theme_classic() +
  labs(x = "Timepoint (HPI)", y = "Protein Abundance", color = "Protein")
  #facet_wrap(~ category, scales = "free_y")

# Plot for Botrydial
botrydial_levels <- paste0("Bcbot", 1:5)

df_filt %>%
  filter(category == "Botrydial") %>%
  filter(gene_name %in% botrydial_levels) %>%
  mutate(gene_name = factor(gene_name, levels = botrydial_levels)) %>%
  ggplot(aes(x = hpi, y = abundance, color = gene_name)) +
  geom_boxplot() +
  scale_color_brewer(palette = "Paired") +
  theme_classic() +
  labs(x = "Timepoint (HPI)", y = "Protein Abundance", color = "Protein")
  #facet_wrap(~ category, scales = "free_y")

# Plot for Polygalacturonase
polygalacturonase_levels <- unique(df_filt %>%
  filter(category == "Polygalacturonase") %>%
  pull(gene_name))

df_filt %>%
  filter(category == "Polygalacturonase") %>%
  filter(gene_name %in% polygalacturonase_levels) %>%
  mutate(gene_name = factor(gene_name, levels = polygalacturonase_levels)) %>%
  ggplot(aes(x = hpi, y = abundance)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Timepoint (HPI)", y = "Protein Abundance")
  #facet_wrap(~ category, scales = "free_y")

