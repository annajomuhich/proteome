### Summarizing total number of significant proteins for each Arabidopsis mutant model
### March 2026 AJM

library(tidyverse)
library(ggVennDiagram)

# #add column for model name
# at_mybsvsCol0_anova <- at_mybsvsCol0_anova %>%
# 	mutate(model = "myb28/29 vs Col0")
# at_aop2vsCol0_anova <- at_aop2vsCol0_anova %>%
# 	mutate(model = "AOP2 vs Col0")
# at_tgg12vsCol0_anova <- at_tgg12vsCol0_anova %>%
# 	mutate(model = "tgg12 vs Col0")
# at_accessions_anova <- at_accessions_anova %>%
# 	mutate(model = "WT accessions")

tgg <- read.csv("data/mutant/mut_at_model_20260330/at_tgg12vsCol0/at_tgg12vsCol0_anova.csv")
aop <- read.csv("data/mutant/mut_at_model_20260330/at_aop2vsCol0/at_aop2vsCol0_anova.csv")
myb <- read.csv("data/mutant/mut_at_model_20260330/at_mybsvsCol0/at_mybsvsCol0_anova.csv")
acc <- read.csv("data/mutant/mut_at_model_20260330/at_accessions/at_accessions_anova.csv")

# Create a list of the dataframes and corresponding model names
model_dfs <- list(
	tgg = tgg,
	aop = aop,
	myb = myb,
	acc = acc
)
# For each dataframe, tabulate convergence_note and add a model identifier
convergence_notes_summary <- map2_dfr(
	model_dfs,
	names(model_dfs),
	~ table(.x$convergence_note) %>%
		as.data.frame() %>%
		mutate(model = .y)
) %>%
	pivot_wider(names_from = model, values_from = Freq)
convergence_notes_summary

# Filter each dataframe in model_dfs to rows where convergence_note == "relative convergence (4)"
model_dfs <- map(model_dfs, ~ .x %>%
								 	filter(variable != "intercept") %>%
								 	filter(convergence_note == "relative convergence (4)") %>%
									dplyr::select(protein_ID, variable, p_adj) %>%
									mutate(variable = recode(variable,
																	 `genotype:treatment` = "Genotype X\nTreatment",
																	 treatment = "Infected",
																	 genotype = "Genotype"))
)


# Show the overlap for each variable within each model (not across models)

# For each model, get all significant protein_IDs for each variable
model_variable_lists <- map2(model_dfs, names(model_dfs), function(df, name) {
	# Get list of variables with at least one significant protein
	vars <- unique(df$variable)
	for (var in vars) {
		# For each variable in this model, get protein_IDs that are significant
		sig_proteins <- df %>%
			filter(variable == var, p_adj < 0.05) %>%
			pull(protein_ID) %>%
			unique()
		
		assign(var, sig_proteins, envir = parent.frame())
	}
	# Build a list of variable: protein_IDs for this model
	result <- split(df, df$variable) %>%
		map(~ .x %>%
					filter(p_adj < 0.05) %>%
					pull(protein_ID) %>%
					unique())
	names(result) <- vars
	list(model = name, variable_lists = result)
})

# Define mapping from model_name to desired plot label
model_labels <- c(
	acc = "Accession",
	tgg = "tgg1/2",
	aop = "AOP2",
	myb = "myb28/29"
)

# Plot a Venn diagram for the variables in each model, showing overlap of significant proteins across variables *within* that model
for (entry in model_variable_lists) {
	model_name <- entry$model
	var_lists <- entry$variable_lists
	
	# Look up appropriate label; fall back to model_name if not present in mapping
	plot_label <- if (model_name %in% names(model_labels)) model_labels[[model_name]] else model_name
	
	# Only plot if more than one variable has significant proteins
	num_nonzero <- sum(lengths(var_lists) > 0)
	if (num_nonzero > 1) {
		p <- ggVennDiagram(var_lists, label = "count") +
			scale_fill_gradient(low = "white", high = "purple") +
			theme_void() +
			ggtitle(paste0(plot_label)) +
			coord_cartesian(clip = "off") +
			theme(
				legend.position = "none",
				text = element_text(size = 10),
				plot.title = element_text(hjust = 0.5, face = "bold"),
				plot.margin = margin(20, 40, 20, 40)
			)
		print(p)
		# Optionally, save each plot
		ggsave(paste0("figures/mutant/arabidopsis/venn/mut_at_venn_", model_name, ".png"),
		       plot = p, height = 4, width = 4)
	}
}

#write the significant_genes to a dataframe
# Create a list of all unique genes across all vectors
all_genes <- unique(unlist(significant_genes))

# Create a data frame where each gene has a column for each vector,
# marked "yes" if present, "no" otherwise
significant_genes_df <- data.frame(gene_ID = all_genes)

for (vec_name in names(significant_genes)) {
	significant_genes_df[[vec_name]] <- ifelse(
		significant_genes_df$gene_ID %in% significant_genes[[vec_name]], "yes", "no"
	)
}
significant_genes_df %>% write.csv("data/mutant/mut_bc_model_20260323/model_sig_2XFC_list.csv", row.names = F)




