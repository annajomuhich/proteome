### Mutant data reformat
### March 2026 AJM 

library(tidyverse)
library(readxl)

### Botrytis -------------------------------------
df_raw <- readxl::read_xlsx("data/mutant/input/Botrytis_combined_protein.tsv.xlsx")

# remove contaminant proteins
df <- df_raw %>%
	filter(Organism == "Botryotinia fuckeliana (strain B05.10)")

#remove extra columns
df <- df %>%
	dplyr::select(!c(...57, ...58, ...59, ...60, ...61, `Indistinguishable Proteins`))

# Remove rows where all columns (except columns 1-8) contain 0
cols_to_check <- setdiff(seq_along(df), 1:8)
df <- df[!apply(df[, cols_to_check], 1, function(x) all(x == 0, na.rm = TRUE)), ]

#pull count of Botrytis proteins
#This is the number of proteins in the raw data
df %>% nrow()

#get sample metadata
df <- df %>%
	pivot_longer(cols = all_of(cols_to_check),
							 names_to = "sample_ID",
							 values_to = "abundance") %>%
	mutate(
		genotype = str_extract(sample_ID, "^[^_]+"),
		treatment = if_else(str_detect(sample_ID, "(?<=_)[B]"), "infected", "mock"),
		rep = str_extract(sample_ID, "\\d+$")
	)

#Make sure everything has a Bcin ID.
#Did this by looking up the Protein IDs on Uniprot:
bc_ids <- df %>% pull(`Protein ID`) %>%
	unique()
bc_ids_df <- as.data.frame(bc_ids)
#bc_ids_df %>% write.csv("data/gene_descriptions/Bc_uniprot_IDs.csv", row.names = F)
uniprot <- read.delim("data/gene_descriptions/Bc_uniprot_IDs_output.tsv")
uniprot <- uniprot %>% dplyr::select(Entry, Gene.Names)
uniprot$Bcin <- str_extract(uniprot$Gene.Names, "BCIN_\\d{2}g\\d{5}")
uniprot <- uniprot %>% dplyr::select(-Gene.Names)
uniprot <- uniprot %>%
	mutate(Bcin = gsub(pattern = "BCIN", replacement = "Bcin", Bcin),
				 Bcin = gsub("_", "", Bcin))
uniprot <- uniprot %>%
	dplyr::rename(gene_ID = Bcin)
df <- left_join(df, uniprot, join_by("Protein ID" == "Entry")) %>%
	dplyr::select(gene_ID, everything())

#5 genes are still missing a Bcin ID:
missing <- df %>% filter(is.na(gene_ID))
#fill in the missing ones manually
#can't find reliable IDs for the ones beginning with R9U but I'll just use the BC ID
df <- df %>%
	mutate(gene_ID = if_else(`Protein ID` == "A0A384JAW2", "Bcin03g00480", gene_ID)) %>%
	mutate(gene_ID = if_else(`Protein ID` == "R9U1R1", Gene, gene_ID)) %>%
	mutate(gene_ID = if_else(`Protein ID` == "R9U348", Gene, gene_ID)) %>%
	mutate(gene_ID = if_else(`Protein ID` == "R9U4F5", Gene, gene_ID)) %>%
	mutate(gene_ID = if_else(`Protein ID` == "R9U4K7", Gene, gene_ID))
any(is.na(df$gene_ID))	

#Reorder columns
df <- df %>%
	dplyr::select(gene_ID, sample_ID, abundance, genotype, treatment, rep, everything())

#arrange by gene ID
df <- df %>%
	arrange(gene_ID)

#Some entries of AOP2 have a typo
df <- df %>%
	mutate(genotype = if_else(genotype == "A0P2", "AOP2", genotype))

#keep only the infected samples
df <- df %>%
	filter(treatment == "infected")

# If 2 bioreps of a given sample are NA but the other one has signal, convert that value to NA
# It is probably stochastic detection and not real biology
df <- df %>%
	group_by(gene_ID, treatment, genotype) %>%
	mutate(
		n_detected = sum(abundance == 0, na.rm = TRUE),
		abundance = ifelse(n_detected == 2, 0, abundance)
	) %>%
	ungroup() %>%
	dplyr::select(-n_detected)

#Because we converted more values to NA, some proteins are now entirely NA:
df %>%
	group_by(gene_ID) %>%
	summarise(all_zero = all(abundance == 0)) %>%
	summarise(n_removed = sum(all_zero))
#Removing these:
df <- df %>%
	group_by(gene_ID) %>%
	filter(any(abundance != 0)) %>%
	ungroup()



### Low protein abundance filtering
#Get a list of proteins to keep. These are proteins that:
#For any combo of protein x treatment, there is signal in >= 2 reps
proteins_keep <- df %>%
	filter(abundance != 0) %>%
	group_by(gene_ID, treatment) %>%
	summarise(n_rep = n_distinct(rep), .groups = "drop") %>%
	filter(n_rep >= 2) %>%
	pull(gene_ID) %>%
	unique()

#remove low abundance proteins
df <- df%>%
	filter(gene_ID %in% proteins_keep)

#this protein is a duplicate, let's remove it
df <- df %>%
	filter(`Protein ID` != "A0A384JXJ9")

#count remaining proteins
df %>% pull(gene_ID) %>% unique() %>% length()

df %>%
	write.csv("data/mutant/input/Bc_Proteome_Mut_filtered.csv", row.names = F)

### ------------ Arabidopsis ------------------------------------------------
df_raw <- readxl::read_xlsx("data/mutant/input/Arabidopsis Proteome Quantification.xlsx")

# remove contaminant proteins
df <- df_raw %>%
	filter(Organism == "Arabidopsis thaliana")

#remove extra columns
df <- df %>%
	dplyr::select(!c(Botrytis, Control, Disease, Col, MYB, GLS, `Indistinguishable Proteins`))

# Remove rows where all columns (except columns 1-8) contain 0
cols_to_check <- setdiff(seq_along(df), 1:8)
df <- df[!apply(df[, cols_to_check], 1, function(x) all(x == 0, na.rm = TRUE)), ]

#pull count of Arabidopsis proteins
#This is the number of proteins in the raw data
df %>% nrow()

#get sample metadata
df <- df %>%
	pivot_longer(cols = all_of(cols_to_check),
							 names_to = "sample_ID",
							 values_to = "abundance") %>%
	mutate(
		genotype = str_extract(sample_ID, "^[^_]+"),
		treatment = if_else(str_detect(sample_ID, "(?<=_)[B]"), "infected", "mock"),
		rep = str_extract(sample_ID, "\\d+$")
	)

#Make sure everything has a AT gene ID.
#Did this by looking up the Protein IDs on Uniprot:
at_ids <- df %>% pull(`Protein ID`) %>%
	unique()
at_ids_df <- as.data.frame(at_ids)
#at_ids_df %>% write.csv("data/gene_descriptions/At_uniprot_IDs.csv", row.names = F)
uniprot <- read.delim("data/gene_descriptions/At_uniprot_IDs_output.tsv")
uniprot <- uniprot %>% dplyr::select(Entry, Gene.Names)
uniprot$Ath <- str_extract(uniprot$Gene.Names, "At[0-9CM]g\\d{5}")
uniprot <- uniprot %>% dplyr::select(-Gene.Names)
uniprot <- uniprot %>%
	dplyr::rename(gene_ID = Ath)
df <- left_join(df, uniprot, join_by("Protein ID" == "Entry")) %>%
	dplyr::select(gene_ID, everything())

#Reorder columns
df <- df %>%
	dplyr::select(gene_ID, sample_ID, abundance, genotype, treatment, rep, everything())

#arrange by gene ID
df <- df %>%
	arrange(gene_ID)

#Some entries of AOP2 have a typo
df <- df %>%
	mutate(genotype = if_else(genotype == "A0P2", "AOP2", genotype))

# If 2 bioreps of a given sample are NA but the other one has signal, convert that value to NA
# It is probably stochastic detection and not real biology
df %>% filter(gene_ID == "At1g01220" & genotype == "Sha")
df <- df %>%
	group_by(gene_ID, treatment, genotype) %>%
	mutate(
		n_detected = sum(abundance == 0, na.rm = TRUE),
		abundance = ifelse(n_detected == 2, 0, abundance)
	) %>%
	ungroup() %>%
	dplyr::select(-n_detected)
df %>% filter(gene_ID == "At1g01220" & genotype == "Sha")

#Because we converted more values to NA, some proteins are now entirely NA:
df %>%
	group_by(gene_ID) %>%
	summarise(all_zero = all(abundance == 0)) %>%
	summarise(n_removed = sum(all_zero))
#Removing these:
df <- df %>%
	group_by(gene_ID) %>%
	filter(any(abundance != 0)) %>%
	ungroup()


### Low protein abundance filtering
#Get a list of proteins to keep. These are proteins that:
#For any combo of protein x treatment, there is signal in >= 2 reps
proteins_keep <- df %>%
	filter(abundance != 0) %>%
	group_by(gene_ID, treatment) %>%
	summarise(n_rep = n_distinct(rep), .groups = "drop") %>%
	filter(n_rep >= 2) %>%
	pull(gene_ID) %>%
	unique()

#remove low abundance proteins
df <- df%>%
	filter(gene_ID %in% proteins_keep)

#count remaining proteins
df %>% pull(gene_ID) %>% unique() %>% length()

#reformat Gene IDs
df <- df %>%
	mutate(gene_ID = gsub(pattern = "t", replacement = "T", x = gene_ID)) %>%
	mutate(gene_ID = gsub(pattern = "g", replacement = "G", x = gene_ID))

#write out
df %>% write.csv("data/mutant/input/At_Proteome_Mut_filtered.csv", row.names = F)