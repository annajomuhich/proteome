### Summarizing total number of significant proteins for each Arabidopsis mutant model
### March 2026 AJM

library(tidyverse)

#load anovas for each mutant model
at_mybsvsCol0_anova <- read.csv("data/mutant/mut_at_model_20260330/at_mybsvsCol0/at_mybsvsCol0_anova.csv")
at_aop2vsCol0_anova <- read.csv("data/mutant/mut_at_model_20260330/at_aop2vsCol0/at_aop2vsCol0_anova.csv")
at_tgg12vsCol0_anova <- read.csv("data/mutant/mut_at_model_20260330/at_tgg12vsCol0/at_tgg12vsCol0_anova.csv")
at_accessions_anova <- read.csv("data/mutant/mut_at_model_20260330/at_accessions/at_accessions_anova.csv")

#add column for model name
at_mybsvsCol0_anova <- at_mybsvsCol0_anova %>%
	mutate(model = "myb28/29 vs Col0")
at_aop2vsCol0_anova <- at_aop2vsCol0_anova %>%
	mutate(model = "AOP2 vs Col0")
at_tgg12vsCol0_anova <- at_tgg12vsCol0_anova %>%
	mutate(model = "tgg12 vs Col0")
at_accessions_anova <- at_accessions_anova %>%
	mutate(model = "WT accessions")

#bind anovas into a single dataframe
at_anova <- bind_rows(at_mybsvsCol0_anova, at_aop2vsCol0_anova, at_tgg12vsCol0_anova, at_accessions_anova) %>%
    filter(variable != "intercept")

#count number of significant proteins at each FDR threshold
at_anova_significant <- at_anova %>%
    group_by(model) %>%
	filter(p_adj < 0.05) %>%
	count() %>%
	ungroup()

#dot plot
at_anova_significant %>%
	mutate(variable = case_when(
		variable == "genotype" ~ "Genotype",
		variable == "treatment" ~ "Infected",
		variable == "genotype:treatment" ~ "Genotype x Infected",
		TRUE ~ variable
	)) %>%
    mutate(variable = factor(variable, levels = c("Genotype", "Infected", "Genotype x Infected")),
        model = factor(model, levels = c("myb28/29 vs Col0", "AOP2 vs Col0", "tgg12 vs Col0", "WT accessions"))) %>%
ggplot(aes(x = variable, y = model)) +
  geom_point(aes(size = n)) +
  scale_size(range = c(3, 10)) +
  #scale_fill_manual(values = c("Genotype" = "#155799", "Treatment" = "#155799", "Genotype x Treatment" = "#155799")) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "",
    y = "",
    size = "Significant\nArabidopsis\nproteins"
  ) +
  scale_y_discrete(limits = rev) 
ggsave("figures/mutant/arabidopsis/mut_at_modeldots.png",
			 height = 3.5, width = 4.25)
	



