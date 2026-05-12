### Plotting protein accumulation of Botrytis network members in the GSL AtBc network
### Infecting different arabidopsis mutants

library(tidyverse)

net <- read.csv("../ITCase/input/network/atbc_d50_72_annotated.csv")
prot <- read.csv("data/mutant/input/Bc_Proteome_Mut_filtered.csv") %>%
	mutate(genotype = factor(genotype, levels = c("Col0", "Ler", "Sha","AOP2", "myb28", "myb29", "myb2829", "tgg12"))) %>%
	filter(gene_ID %in% genes$gene) %>%
	dplyr::select(gene_ID, abundance, genotype, rep, Gene, Description)

genes <- net %>%
	filter(
		!startsWith(gene, "AT") &
		gene != "Bcin06g00024" &
		!startsWith(gene_name, "Bcbot")
		) %>%
		dplyr::select(gene, gene_name, gene_desc)

df <- left_join(prot, genes, join_by(gene_ID == gene))

df <- df %>%
	mutate(label = if_else(gene_desc == "N/A", NA,
	if_else(grepl("Zinc", gene_desc), "Alcohol dehydrogenase",
	gene_desc)))

df %>%
	ggplot(aes(x = gene_ID, y = abundance, fill = genotype)) +
	geom_boxplot(size = 0.2) +   # reduced box outline thickness
	theme_classic() +
	labs(x = "", y = "Raw Protein Abundance") +
	scale_fill_brewer(palette = "Accent") +
	scale_x_discrete(
		labels = function(x) {
			custom_labels <- df %>%
				distinct(gene_ID, label) %>%
				filter(gene_ID %in% x)
			setNames(
				paste0(custom_labels$gene_ID, "\n", ifelse(is.na(custom_labels$label), "", custom_labels$label)),
				custom_labels$gene_ID
			)[x]
		}
	) +
	theme(
		axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
	)

ggsave("figures/mutant/botrytis/mut_bc_netmembers.png", height = 5, width = 6,
			 dpi = 2500)
