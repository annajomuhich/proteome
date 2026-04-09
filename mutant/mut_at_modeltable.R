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
    group_by(model, variable) %>%
	filter(p_adj < 0.05) %>%
	count() %>%
	ungroup() %>%
	mutate(variable = recode(variable,
													 `genotype:treatment` = "Genotype X Treatment",
													 treatment = "Infected",
													 genotype = "Genotype")) %>%
	pivot_wider(names_from = variable,
							values_from = n) %>%
	mutate(Total = rowSums(select(., -model))) %>%
	arrange(desc(Total))
at_anova_significant %>% write.csv("tables/mut_at_modeltable.csv", row.names = F)

	



