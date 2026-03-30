### Plotting levels of specific Bcin proteins
### AJM March 2026

library(tidyverse)

anova <- read.csv("tables/mut_bc_modeltable.csv")
expr <- read.csv("data/mutant/input/Bc_Proteome_Mut_filtered.csv")

### ----------------- 5 all-mutant protein levels -------------------------

genes <- anova %>%
	filter(tgg < 0.05 &
				 	aop < 0.05 &
				 	myb < 0.05 &
				 	acc < 0.05
	) %>%
	pull(protein_ID)

genes5 <- expr %>%
  filter(gene_ID %in% genes & treatment == "infected") %>%
  group_by(gene_ID, genotype) %>%
  summarise(
    abundance_mean = mean(abundance),
    abundance_sd = sd(abundance),
    .groups = "drop"
  ) %>%
  mutate(genotype = factor(genotype, levels = c("Col0", "Ler", "Sha","AOP2", "myb28", "myb29", "myb2829", "tgg12")))

genes5 %>%
	mutate(genotype = dplyr::recode(genotype,
																	"tgg12" = "tgg1/2",
																	"myb2829" = "myb28/29"
	)) %>%
	mutate(gene_desc = if_else(gene_ID == "Bcin02g06940", "Bcgas1",
		if_else(gene_ID == "Bcin06g00024", "Beta-lactamase",
		if_else(gene_ID == "Bcin06g00620", "Bctpp2",
		if_else(gene_ID == "Bcin11g06240", "Bcpgk1",
		if_else(gene_ID == "Bcin13g05810", "Aldehyde dehydrogenase", "N/A")))))) %>%
	mutate(label = paste0(gene_ID, "\n", gene_desc)) %>%
	ggplot(aes(x = genotype, y = abundance_mean, fill = genotype)) +
	geom_col(width = 0.7, show.legend = FALSE) +
	geom_errorbar(
		aes(ymin = abundance_mean - abundance_sd, ymax = abundance_mean + abundance_sd),
		width = 0.2
	) +
	labs(x = "", y = "Mean Abundance") +
	facet_wrap(~ label, scales = "free_y") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_brewer(palette = "Accent")
ggsave("figures/mutant/botrytis/mut_bc_5allmutant_proteinlevels.png",
			 height = 4.5, width = 7)

### ------------- plots for significant >2XFC proteins ---------------

df <- read.csv("data/mutant/mut_bc_model_20260323/model_sig_2XFC_list.csv")
annot <- read.csv("data/gene_descriptions/Bcin_Annotations_Full_transcript.csv") %>%
	dplyr::select(X.Gene.ID., X.Gene.Name.or.Symbol., X.PFam.Description.) %>%
	dplyr::rename(protein_ID = X.Gene.ID.,
								gene_name = X.Gene.Name.or.Symbol.,
								gene_desc = X.PFam.Description.)

#Any overlap plot
overlap <- df %>%
	filter(rowSums(across(where(is.character), ~ .x == "yes")) >= 2) %>%
	pull(gene_ID)
overlap_annot <- annot %>% filter(protein_ID %in% overlap)
overlap_expr <- expr %>%
	filter(gene_ID %in% overlap & treatment == "infected") %>%
	group_by(gene_ID, genotype) %>%
	summarise(
		abundance_mean = mean(abundance),
		abundance_sd = sd(abundance),
		.groups = "drop"
	) %>%
	mutate(genotype = factor(genotype, levels = c("Col0", "Ler", "Sha","AOP2", "myb28", "myb29", "myb2829", "tgg12")))

overlap_expr %>%
	mutate(genotype = dplyr::recode(genotype,
																	"tgg12" = "tgg1/2",
																	"myb2829" = "myb28/29"
	)) %>%
	mutate(gene_desc = if_else(gene_ID == "Bcin08g00280", "Serine carboxypeptidase",
														 if_else(gene_ID == "Bcin06g00024", "Beta-lactamase",
														 				if_else(gene_ID == "Bcin10g03910", "Bchsp104",
														 								if_else(gene_ID == "Bcin13g00010", "CYP450",
														 												if_else(gene_ID == "Bcin15g04630", "Bcnsr1", "N/A")))))) %>%
	mutate(label = paste0(gene_ID, "\n", gene_desc)) %>%
	ggplot(aes(x = genotype, y = abundance_mean, fill = genotype)) +
	geom_col(width = 0.7, show.legend = FALSE) +
	geom_errorbar(
		aes(ymin = abundance_mean - abundance_sd, ymax = abundance_mean + abundance_sd),
		width = 0.2
	) +
	labs(x = "", y = "Mean Abundance") +
	facet_wrap(~ label, scales = "free_y") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_brewer(palette = "Accent")
ggsave("figures/mutant/botrytis/mut_bc_5overlap2XFC_proteinlevels.png",
			 height = 4.5, width = 7)

#tgg1/2 only plot
tgg12 <- df %>%
	filter(tgg1.2 == "yes" &
				 aop2 == "no" &
				 myb28.29 == "no" &
				 Accession == "no") %>%
	pull(gene_ID)
tgg12_expr <- expr %>%
	filter(gene_ID %in% tgg12 & treatment == "infected") %>%
	group_by(gene_ID, genotype) %>%
	summarise(
		abundance_mean = mean(abundance),
		abundance_sd = sd(abundance),
		.groups = "drop"
	) %>%
	filter(genotype == "Col0" | genotype == "tgg12") %>%
	mutate(genotype = factor(genotype, 
													 levels = c("Col0", "tgg12")))
tgg12_expr %>%
	mutate(genotype = dplyr::recode(genotype,
																	"tgg12" = "tgg1/2"
	)) %>%
	mutate(gene_desc = if_else(gene_ID == "Bcin13g05240", "Bcyra1", "N/A")) %>%
	mutate(label = paste0(gene_ID, "\n", gene_desc)) %>%
	ggplot(aes(x = genotype, y = abundance_mean, fill = genotype)) +
	geom_col(width = 0.7, show.legend = FALSE) +
	geom_errorbar(
		aes(ymin = abundance_mean - abundance_sd, ymax = abundance_mean + abundance_sd),
		width = 0.2
	) +
	labs(x = "", y = "Mean Abundance") +
	facet_wrap(~ label, scales = "free_y") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(
		values = scales::brewer_pal(palette = "Accent")(8)[c(1, 8)],
		limits = c("Col0", "tgg1/2")
	)
ggsave("figures/mutant/botrytis/mut_bc_tgg12only_proteinlevels.png",
			 height = 3, width = 1.7, dpi = 2000)

#AOP2 only plot
aop2 <- df %>%
	filter(aop2 == "yes" &
				 tgg1.2 == "no" &
				 myb28.29 == "no" &
				 Accession == "no") %>%
	pull(gene_ID)
aop2_annot <- annot %>%
	filter(protein_ID %in% aop2)
aop2_expr <- expr %>%
	filter(gene_ID %in% aop2 & treatment == "infected") %>%
	group_by(gene_ID, genotype) %>%
	summarise(
		abundance_mean = mean(abundance),
		abundance_sd = sd(abundance),
		.groups = "drop"
	) %>%
	filter(genotype == "Col0" | genotype == "AOP2") %>%
	mutate(genotype = factor(genotype, 
													 levels = c("Col0", "AOP2")))
aop2_expr %>%
	mutate(genotype = dplyr::recode(genotype,
																	"AOP2" = "AOP2"
	)) %>%
	mutate(gene_desc = if_else(gene_ID == "Bcin08g02080", "Antibiotic\nbiosynthesis\nmonooxygenase", 
														 if_else(gene_ID == "Bcin11g05080", "Aldo/keto\nreductase", 
														 				if_else(gene_ID == "Bcin14g00610", "Bcpg2",
														 								if_else(gene_ID == "Bcin14g05500", "Bcgod1", "N/A"))))) %>%
	mutate(label = paste0(gene_ID, "\n", gene_desc)) %>%
	ggplot(aes(x = genotype, y = abundance_mean, fill = genotype)) +
	geom_col(width = 0.7, show.legend = FALSE) +
	geom_errorbar(
		aes(ymin = abundance_mean - abundance_sd, ymax = abundance_mean + abundance_sd),
		width = 0.2
	) +
	labs(x = "", y = "Mean Abundance") +
	facet_wrap(~ label, scales = "free_y", nrow = 1) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		scale_fill_manual(
		values = scales::brewer_pal(palette = "Accent")(8)[c(1, 4)],
		limits = c("Col0", "AOP2")
	)
ggsave("figures/mutant/botrytis/mut_bc_aop2only_proteinlevels.png",
			 height = 3, width = 6, dpi = 2000)

#myb28/29 only plot
myb2829 <- df %>%
	filter(myb28.29 == "yes" &
				 tgg1.2 == "no" &
				 aop2 == "no" &
				 Accession == "no") %>%
	pull(gene_ID)
myb2829_annot <- annot %>%
	filter(protein_ID %in% myb2829)
myb2829_expr <- expr %>%
	filter(gene_ID %in% myb2829 & treatment == "infected") %>%
	group_by(gene_ID, genotype) %>%
	summarise(
		abundance_mean = mean(abundance),
		abundance_sd = sd(abundance),
		.groups = "drop"
	) %>%
	filter(genotype == "Col0" | genotype == "myb28" | genotype == "myb29" | genotype == "myb2829") %>%
	mutate(genotype = factor(genotype, 
													 levels = c("Col0", "myb28", "myb29", "myb2829")))
myb2829_expr <- left_join(myb2829_expr,
													read.csv("data/genes_of_interest/mut_bc_myb_sig2XFC_genes_20260325.csv"),
													by = "gene_ID")
myb2829_expr %>%
	mutate(genotype = dplyr::recode(genotype,
																	"myb2829" = "myb28/29"
	)) %>%
	mutate(gene_desc = gsub("\\\\n", "\n", gene_desc)) %>%
	mutate(label = paste0(gene_ID, "\n", gene_desc)) %>%
	ggplot(aes(x = genotype, y = abundance_mean, fill = genotype)) +
	geom_col(width = 0.7, show.legend = FALSE) +
	geom_errorbar(
		aes(ymin = abundance_mean - abundance_sd, ymax = abundance_mean + abundance_sd),
		width = 0.2
	) +
	labs(x = "", y = "Mean Abundance") +
	facet_wrap(~ label, scales = "free_y", nrow = 4) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(
		values = c("#7fc97f", "#386cb0", "#f0017f", "#bf5b17"),
		limits = c("Col0", "myb28", "myb29", "myb28/29")
	)
ggsave("figures/mutant/botrytis/mut_bc_mybonly_proteinlevels.png",
			 height = 7, width = 3, dpi = 2000)

#Accession only plot
accession <- df %>%
	filter(Accession == "yes" &
				 tgg1.2 == "no" &
				 aop2 == "no" &
				 myb28.29 == "no") %>%
	pull(gene_ID)
accession_expr <- expr %>%
	filter(gene_ID %in% accession & treatment == "infected") %>%
	group_by(gene_ID, genotype) %>%
	summarise(
		abundance_mean = mean(abundance),
		abundance_sd = sd(abundance),
		.groups = "drop"
	) %>%
	filter(genotype == "Col0" | genotype == "Ler" | genotype == "Sha") %>%	
	mutate(genotype = factor(genotype, 
													 levels = c("Col0", "Ler", "Sha")))
accession_expr <- left_join(accession_expr,
														read.csv("data/genes_of_interest/mut_bc_acc_sig2XFC_genes_20260326.csv"),
														by = "gene_ID")
accession_expr %>%
	mutate(gene_desc = if_else(is.na(gene_desc), "", gene_desc)) %>%
	mutate(gene_desc = gsub("\\\\n", "\n", gene_desc)) %>%
	mutate(label = paste0(gene_ID, "\n", gene_desc)) %>%
	ggplot(aes(x = genotype, y = abundance_mean, fill = genotype)) +
	geom_col(width = 0.7, show.legend = FALSE) +
	geom_errorbar(
		aes(ymin = abundance_mean - abundance_sd, ymax = abundance_mean + abundance_sd),
		width = 0.2
	) +
	labs(x = "", y = "Mean Abundance") +
	facet_wrap(~ label, scales = "free_y") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_manual(
		values = scales::brewer_pal(palette = "Accent")(8)[c(1, 2, 3)],
		limits = c("Col0", "Ler", "Sha")
	)
ggsave("figures/mutant/botrytis/mut_bc_accessiononly_proteinlevels.png",
			 height = 7.5, width = 6.5, dpi = 2000)
