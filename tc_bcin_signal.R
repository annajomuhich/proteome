### Bcin Protein Signal over Time
### February 2026 AJM

library(tidyverse)

df <- read.csv("data/timecourse/input/AtBc_Proteome_TimeCourse_filtered.csv")

samplekey <- df %>%
	select(sample_ID, genotype, treatment, hpi, rep) %>%
	distinct(.keep_all = T)

summary <- df %>%
	filter(treatment == "infected") %>%
	group_by(sample_ID, organism) %>%
	summarise(sum = sum(abundance)) %>%
	pivot_wider(names_from = organism,
							values_from = sum) %>%
	mutate(bcin_pct_abundance = botrytis / (arabidopsis+botrytis) * 100) #%>%

summary <- left_join(summary, samplekey, by = "sample_ID")

summary$hpi <- as.factor(summary$hpi)

summary %>%
	ggplot(aes(x = hpi, y = bcin_pct_abundance, color = genotype)) +
	geom_jitter(size = 2.5, width = 0.15) +
	theme_minimal() +
	ylab("Botrytis protein abundance %") +
	xlab("HPI") +
	scale_color_manual(values = c("red", "blue"))
ggsave("figures/timecourse/botrytis/tc_bcin_signal_pct.png", height = 4, width = 5)

summary %>%
	ggplot(aes(x = hpi, y = botrytis, color = genotype)) +
	geom_jitter(size = 2.5, width = 0.15) +
	theme_minimal() +
	ylab("Botrytis protein abundance") +
	xlab("HPI") +
	scale_color_manual(values = c("red", "blue"))
ggsave("figures/timecourse/botrytis/tc_bcin_signal.png", height = 4, width = 5)
