### Get XP ID key for Botrytis proteins
### March 2026 AJM

gff <- read_tsv("data/gene_descriptions/ncbi_dataset/ncbi_dataset/data/GCF_000143535.2/genomic.gff", comment = "#", col_names = FALSE)

parsed <- gff %>%
	mutate(
		XP = str_extract(X9, "(?<=protein_id=)[^;]+"),
		BCIN = str_extract(X9, "(?<=locus_tag=)[^;]+"),
		gene_name = str_extract(X9, "(?<=gene=)[^;]+"),
		gene_ID = str_replace(BCIN, "^BCIN_", "Bcin")
	)

df <- parsed %>%
	select(XP, gene_ID, gene_name) %>%
	filter(startsWith(XP, "XP")) %>%
	unique()

df %>%
	write.csv("data/gene_descriptions/Bcin_XP_key.csv", row.names = F)

