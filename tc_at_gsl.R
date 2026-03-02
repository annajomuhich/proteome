### At Time Course - GSL
### February 2026 AJM

library(tidyverse)

gsl_df <- read.csv("data/genes_of_interest/gsl_genes.csv")

#define list of gsl genes of interest
glucosinolate_top <- c(
	"MAM1", "MAM3",
	"BCAT4",
	"CYP79F1",
	"CYP83A1",
	"UGT74B"
)
gsl_df <- gsl_df %>% filter(gene_name %in% glucosinolate_top)
gsl<- gsl_df %>% pull(gene)

df <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_adjusted_emmeans_log.csv")

#remove `.1` from proten/gene ID, e.g.
df <- df %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))

#filter to cam genes
df <- df %>%
	filter(protein_ID %in% gsl)

df <- left_join(df, gsl_df, join_by(protein_ID == gene))

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
#		legend.position = NA,   # x, y
		plot.title = element_text(hjust = 0.5, face = "bold")
		#legend.key.size = unit(0.4, "cm")
	) +
	ggtitle("Glucosinolate Proteins") +
	scale_color_brewer(palette = "Dark2", direction = -1)

ggsave("figures/timecourse/arabidopsis/metabolites/gsl.png",
			 height = 3.5, width = 6.5)
