### Visualizing overlap of Botrytis proteins significantly different across mutant comparisons
### March 2026 AJM

library(tidyverse)
library(ggVennDiagram)
# library(eulerr)

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
								 	filter(variable == "genotype") %>%
								 	filter(convergence_note == "relative convergence (4)") 
								 	)

#for each dataframe, get a vector of significant genes
significant_genes <- map(model_dfs, ~ .x %>%
								 	filter(p_adj < 0.05) %>%
								 	pull(protein_ID) %>%
								 	unique()
								 	)

names(significant_genes) <- c("tgg1/2", "AOP2", "myb28/29", "Accession")

### Plot all overlaps ----------------------------------------------------------
ggVennDiagram(significant_genes, label = "count") +
	scale_fill_gradient(low = "white", high = "purple") +  # controls region shading
	# scale_color_manual(values = cols) +                       # set outlines
	theme_void() +
	coord_cartesian(clip = "off") +   #allow drawing outside panel
	theme(
		legend.position = "none",
		text = element_text(size = 10),
			plot.margin = margin(20, 40, 20, 40)  #add space (top, right, bottom, left)
	)
#ggsave("figures/mutant/botrytis/mut_bc_venn.png", height = 4, width = 4)

### Plot tgg/mybs and aop/acc separately -------------------------------------
tggmyb <- list(significant_genes[["tgg1/2"]], significant_genes[["myb28/29"]])
names(tggmyb) <- c("tgg1/2", "myb28/29")
ggVennDiagram(tggmyb, label = "count", label_size = 4.5) +
	scale_fill_gradient(low = "white", high = "grey") +  # controls region shading
	theme_void() +
	coord_flip() +
	#coord_fixed(clip = "off") +   #allow drawing outside panel
	theme(
		legend.position = "none",
		text = element_text(size = 10),
		plot.margin = margin(20, 30, 20, 30)  #add space (top, right, bottom, left)
	)
#ggsave("figures/mutant/botrytis/mut_bc_venn_tggmyb.png", width = 4.12, height = 3)

aopacc <- list(significant_genes[["aop2"]], significant_genes[["Accession"]])
names(aopacc) <- c("AOP2", "Accession")
ggVennDiagram(aopacc, label = "count", label_size = 4.5) +
	scale_fill_gradient(low = "white", high = "grey") +  # controls region shading
	theme_void() +
	coord_flip() +
	#coord_fixed(clip = "off") +   #allow drawing outside panel
	theme(
		legend.position = "none",
		text = element_text(size = 10),
		plot.margin = margin(20, 30, 20, 30)  #add space (top, right, bottom, left)
	)
#ggsave("figures/mutant/botrytis/mut_bc_venn_aopacc.png", width = 4.12, height = 3)

### ----- Trying with a 2X FC filter -------------------------------

significant_genes

deg <- list(
	tgg = read.csv("data/mutant/mut_bc_model_20260323/bc_tgg12vsCol0/bc_tgg12vsCol0_DEGs.csv"),
	aop = read.csv("data/mutant/mut_bc_model_20260323/bc_aop2vsCol0/bc_aop2vsCol0_DEGs.csv"),
	myb = read.csv("data/mutant/mut_bc_model_20260323/bc_mybsvsCol0/bc_mybsvsCol0_DEGs.csv"),
	acc = read.csv("data/mutant/mut_bc_model_20260323/bc_accessions/bc_accessions_DEGs.csv")
)

#filter each dataframe to those with a 2X FC from Col0
deg <- map(deg, ~ .x %>%
	filter(log2FC > 1 | log2FC < -1)
)

#pull out the protein_IDs for each dataframe
deg_ids <- map(deg, ~ .x %>%
	pull(protein_ID)
)

# For each vector in significant_genes, keep only elements also found in the corresponding deg_ids vector
significant_genes <- map2(significant_genes, deg_ids, function(sig, deg) {
	sig[sig %in% deg]
})

#plot the venn diagram
ggVennDiagram(significant_genes, label = "count") +
	scale_fill_gradient(low = "white", high = "purple") +
	theme_void() +
	coord_cartesian(clip = "off") +
	theme(legend.position = "none", text = element_text(size = 10), plot.margin = margin(20, 40, 20, 40))
ggsave("figures/mutant/botrytis/mut_bc_venn_2XFC.png", width = 4, height = 4)

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
