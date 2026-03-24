### Histograms of protein abundance
### March 2026 AJM

df <- read.csv("data/mutant/input/Bc_Proteome_Mut_filtered.csv") %>%
	filter(treatment == "infected")

library(ggplot2)

# # Histogram for arabidopsis
# df_arabidopsis <- subset(df, organism == "arabidopsis" & 
# 												 	abundance >= 0.8 &
# 												 	abundance < 50000)
# ggplot(df_arabidopsis, aes(x = abundance)) +
#   geom_histogram(binwidth = 100, fill = "forestgreen") +
#   theme_minimal() +
#   xlab("Protein Abundance") +
#   ylab("Number of Arabidopsis Proteins")
# ggsave("figures/timecourse/arabidopsis/tc_at_histogram.png", height = 4, width = 4)

# # Histogram for botrytis
df_botrytis <- subset(df, abundance >= 1352065 &
												abundance < 500000000)
ggplot(df_botrytis, aes(x = abundance)) +
  geom_histogram(binwidth = 100000, fill = "orange") +
  theme_minimal() +
  xlab("Protein Abundance") +
  ylab("Number of Botrytis Proteins")
ggsave("figures/mutant/botrytis/mut_bc_histogram.png", height = 4, width = 4)

