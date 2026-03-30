### Making summary table of model results - Botrytis mutant proteome
### March 2026 AJM

library(tidyverse)

tgg <- read.csv("data/mutant/mut_bc_model_20260323/bc_tgg12vsCol0/bc_tgg12vsCol0_anova.csv")
aop <- read.csv("data/mutant/mut_bc_model_20260323/bc_aop2vsCol0/bc_aop2vsCol0_anova.csv")
myb <- read.csv("data/mutant/mut_bc_model_20260323/bc_mybsvsCol0/bc_mybsvsCol0_anova.csv")
acc <- read.csv("data/mutant/mut_bc_model_20260323/bc_accessions/bc_accessions_anova.csv")

# Create a list of the dataframes and corresponding model names
model_dfs <- list(
  tgg = tgg,
  aop = aop,
  myb = myb,
  acc = acc
)

# Filter each dataframe in model_dfs to rows where convergence_note == "relative convergence (4)"
model_dfs <- map(model_dfs, ~ .x %>%
								 	filter(variable == "genotype") %>%
								 	filter(convergence_note == "relative convergence (4)") 
								 	)

#make a combined dataframe of all the model dataframes
df_all <- bind_rows(
  lapply(names(model_dfs), function(name) {
    df <- model_dfs[[name]]
    df$model_name <- name
    df
  })
) %>%
	dplyr::select(protein_ID, p_adj, model_name) %>%
	pivot_wider(names_from = model_name,
							values_from = p_adj)

df_sig <- df_all %>%
  filter(if_any(-protein_ID, ~ .x < 0.05))

annot <- read.csv("data/gene_descriptions/Bcin_Annotations_Full_transcript.csv")
annot <- annot %>%
	dplyr::select(X.Gene.ID., X.Gene.Name.or.Symbol., X.PFam.Description.) %>%
	dplyr::rename(protein_ID = X.Gene.ID.,
								gene_name = X.Gene.Name.or.Symbol.,
								gene_desc = X.PFam.Description.)
df_sig <- left_join(df_sig, annot, by = "protein_ID")
df_sig <- distinct(df_sig) #there's some duplicates now not sure why

#one additional duplciated protein from the annotation. removing
df_sig %>% filter(duplicated(df_sig$protein_ID))
df_sig <- df_sig %>%
  filter(!(protein_ID == "Bcin03g02520" & gene_desc == "N/A"))

df_sig %>% write.csv("tables/mut_bc_modeltable.csv", row.names = F)

df_sig %>%
	filter(tgg < 0.05 &
				 	aop < 0.05 &
				 	myb < 0.05 &
				 	acc < 0.05
	)
