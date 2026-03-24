### At Time Course - GSL
### February 2026 AJM

library(tidyverse)

gsl_df <- read.csv("data/genes_of_interest/gsl_cam_genes.csv")

gsl<- gsl_df %>% filter(Pathways == "Indolic Glucosinolate") %>% pull(AGI)

df <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_adjusted_emmeans_response.csv")

#remove `.1` from proten/gene ID, e.g.
df <- df %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))

#filter to cam genes
df <- df %>%
	filter(protein_ID %in% gsl)

df <- left_join(df, gsl_df, join_by(protein_ID == AGI))

df_centered <- df %>%
	group_by(protein_ID) %>%
	mutate(emmean_centered = emmean_response - mean(emmean_response, na.rm = TRUE)) %>%
	ungroup()

df_centered %>%
	ggplot(aes(x = hpi,
						 y = emmean_centered,
						 color = treatment,
						 group = treatment)) +
	geom_line(size = 1) +
	geom_point(size = 2) +
	facet_wrap(~ GeneName, scales = "free_y") +
	theme_minimal() +
	labs(
		x = "Time (hpi)",
		y = "Protein Abundance",
		color = "Treatment"
	) +
	theme(
#		legend.position = NA,   # x, y
		plot.title = element_text(hjust = 0.5, face = "bold")
		#legend.key.size = unit(0.4, "cm")
	) +
	ggtitle("Indolic Glucosinolate Biosynthetic Proteins") +
	scale_color_manual(values = c("black", "grey")) +
	scale_y_continuous(labels = scales::label_scientific())

ggsave("figures/timecourse/arabidopsis/metabolites/gsl_indolic.png",
			 height = 5.5, width = 7.5)

### aliphatic GSL --------------------------------------------------------
gsl_df <- read.csv("data/genes_of_interest/gsl_cam_genes.csv")

gsl<- gsl_df %>% filter(Pathways == "Aliphatic Glucosinolate") %>% pull(AGI)

df <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_adjusted_emmeans_response.csv")

#remove `.1` from proten/gene ID, e.g.
df <- df %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))

#filter to cam genes
df <- df %>%
	filter(protein_ID %in% gsl)

df <- left_join(df, gsl_df, join_by(protein_ID == AGI))

df_centered <- df %>%
	group_by(protein_ID) %>%
	mutate(emmean_centered = emmean_response - mean(emmean_response, na.rm = TRUE)) %>%
	ungroup()

df_centered %>%
	ggplot(aes(x = hpi,
						 y = emmean_centered,
						 color = treatment,
						 group = treatment)) +
	geom_line(size = 1) +
	geom_point(size = 2) +
	facet_wrap(~ GeneName, scales = "free_y") +
	theme_minimal() +
	labs(
		x = "Time (hpi)",
		y = "Protein Abundance",
		color = "Treatment"
	) +
	theme(
		#		legend.position = NA,   # x, y
		plot.title = element_text(hjust = 0.5, face = "bold")
		#legend.key.size = unit(0.4, "cm")
	) +
	ggtitle("Aliphatic Glucosinolate Biosynthetic Proteins") +
	scale_color_manual(values = c("black", "grey")) +
	scale_y_continuous(labels = scales::label_scientific())

ggsave("figures/timecourse/arabidopsis/metabolites/gsl_aliphatic.png",
			 height = 8.5, width = 10.5)

### GSL activation --------------------------------------------------------
gsl_df <- read.csv("data/genes_of_interest/gsl_cam_genes.csv")

gsl<- gsl_df %>% filter(Pathways == "Activation") %>% pull(AGI)

df <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_adjusted_emmeans_response.csv")

#remove `.1` from proten/gene ID, e.g.
df <- df %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))

#filter to cam genes
df <- df %>%
	filter(protein_ID %in% gsl)

df <- left_join(df, gsl_df, join_by(protein_ID == AGI))

df_centered <- df %>%
	group_by(protein_ID) %>%
	mutate(emmean_centered = emmean_response - mean(emmean_response, na.rm = TRUE)) %>%
	ungroup()

df_centered %>%
	ggplot(aes(x = hpi,
						 y = emmean_centered,
						 color = treatment,
						 group = treatment)) +
	geom_line(size = 1) +
	geom_point(size = 2) +
	facet_wrap(~ GeneName, scales = "free_y") +
	theme_minimal() +
	labs(
		x = "Time (hpi)",
		y = "Protein Abundance",
		color = "Treatment"
	) +
	theme(
		#		legend.position = NA,   # x, y
		plot.title = element_text(hjust = 0.5, face = "bold")
		#legend.key.size = unit(0.4, "cm")
	) +
	ggtitle("GSL Activation Proteins") +
	scale_color_manual(values = c("black", "grey")) +
	scale_y_continuous(labels = scales::label_scientific())

ggsave("figures/timecourse/arabidopsis/metabolites/gsl_activation.png",
			 height = 5.5, width = 7.5)

### Beta-glucosidases --------------------------------------------------------
gsl_df <- read.csv("data/genes_of_interest/bglucosidases.csv")

gsl<- gsl_df  %>% pull(AGI)

df <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_adjusted_emmeans_response.csv")

#remove `.1` from proten/gene ID, e.g.
df <- df %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))

#filter to cam genes
df <- df %>%
	filter(protein_ID %in% gsl)

df <- left_join(df, gsl_df, join_by(protein_ID == AGI))

df_centered <- df %>%
	group_by(protein_ID) %>%
	mutate(emmean_centered = emmean_response - mean(emmean_response, na.rm = TRUE)) %>%
	ungroup()

df_centered %>%
	ggplot(aes(x = hpi,
						 y = emmean_centered,
						 color = treatment,
						 group = treatment)) +
	geom_line(size = 1) +
	geom_point(size = 2) +
	facet_wrap(~ GeneName, scales = "free_y") +
	theme_minimal() +
	labs(
		x = "Time (hpi)",
		y = "Protein Abundance",
		color = "Treatment"
	) +
	theme(
		#		legend.position = NA,   # x, y
		plot.title = element_text(hjust = 0.5, face = "bold")
		#legend.key.size = unit(0.4, "cm")
	) +
	ggtitle("Beta-glucosidases") +
	scale_color_manual(values = c("black", "grey")) +
	scale_y_continuous(labels = scales::label_scientific())

ggsave("figures/timecourse/arabidopsis/metabolites/bgals.png",
			 height = 5.5, width = 7.5)

