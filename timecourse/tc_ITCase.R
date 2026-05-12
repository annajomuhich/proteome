### Plotting Bcin06g00024 abundance over time
### April 2026 AJM

library(tidyverse)

df <- read.csv("data/timecourse/input/AtBc_Proteome_TimeCourse_filtered.csv")

df_filt <- df %>%
    filter(treatment == "infected") %>%
    filter(protein_ID == "Bcin06g00024") %>%
    filter(genotype == "Col.0") %>%
    mutate(hpi = as.factor(hpi)) 

ggplot(df_filt, aes(x = hpi, y = abundance)) +
    geom_boxplot() +
		#geom_jitter() + 
    #geom_line() +
    theme_classic() +
    labs(x = "Timepoint (HPI)", y = "BcSaxA Protein Abundance")

ggsave("figures/timecourse/botrytis/tc_ITCase.png", width = 3, height = 2.5, dpi = 2000)

### With mock

df_filt <- df %>%
	filter(protein_ID == "Bcin06g00024") %>%
	filter(genotype == "Col.0") %>%
	mutate(hpi = as.factor(hpi)) 

ggplot(df_filt, aes(x = hpi, y = abundance, fill = treatment)) +
	geom_boxplot() +
	#geom_jitter() + 
	#geom_line() +
	theme_classic() +
	labs(x = "Timepoint (HPI)", y = "BcSaxA Protein Abundance")
ggsave("figures/timecourse/botrytis/tc_ITCase_wMock.png", width = 3.5, height = 2.5, dpi = 2000)
