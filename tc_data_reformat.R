### Timecourse data reformat
### February 2026 AJM 

library(tidyverse)

df_raw <- read.delim("data/timecourse/input/0466_DIA-FA-GL_Botrytis-timecourse_54Samples_IntMatrix.tsv")

# Remove rows that are entirely NA except for column PG.ProteinGroups
df <- df_raw[!apply(df_raw[, setdiff(names(df_raw), "PG.ProteinGroups")], 1, function(x) all(is.na(x))), ]

# remove contaminant proteins
df <- df %>%
	filter(grepl("^AT", PG.ProteinGroups) |
				 	grepl("^XP", PG.ProteinGroups))

#pull count of Arabidopsis and Botrytis proteins
#This is the number of proteins for each organism in the raw data
df %>%
	mutate(prefix = str_sub(PG.ProteinGroups, 1, 2)) %>%
	count(prefix)

#get sample metadata
df <- df %>%
	pivot_longer(cols = !PG.ProteinGroups,
							 names_to = "sample_ID",
							 values_to = "abundance") %>%
	mutate(
		genotype = str_extract(sample_ID, "^[^_]+"),
		treatment = if_else(str_detect(sample_ID, "(?<=_)[M]"), "mock", "infected"),
		hpi = str_extract(sample_ID, "(?<=_[A-Z])\\d{2}") %>% as.numeric(),
		rep = str_extract(sample_ID, "(?<=_)\\d+$") %>% as.numeric(),
		protein_ID = str_extract(`PG.ProteinGroups`, "^[^;]+"),
		organism = if_else(grepl("^AT", PG.ProteinGroups), "arabidopsis", "botrytis")
	)

#I believe the proteins are isoform groups, grabbed the first one in `protein_ID`
#In case I want a single value

#Reorder columns
df <- df %>%
	dplyr::select(protein_ID, PG.ProteinGroups, organism, sample_ID, genotype, treatment, hpi, rep, abundance)

# #Replace NaN with 0
# df <- df %>%
# 	mutate(abundance = if_else(is.na(abundance), 0, abundance))
#Replace NaN with the lowest observed value
min <- df %>%
	filter(!is.na(abundance)) %>%
	pull(abundance) %>%
	min()
df <- df %>%
	mutate(abundance = if_else(is.na(abundance), min, abundance))

### Low protein abundance filtering
#Get a list of proteins to keep. These are proteins that:
#For any combo of protein x treatment x hpi, there is signal in >= 2 reps
proteins_keep <- df %>%
	filter(abundance > min) %>% #using the min value (~0.71) as threshold
	group_by(protein_ID, treatment, hpi) %>%
	summarise(n_rep = n_distinct(rep), .groups = "drop") %>%
	filter(n_rep >= 2) %>%
	pull(protein_ID) %>%
	unique()

#remove low abundance proteins
df <- df%>%
	filter(protein_ID %in% proteins_keep)

#pull count of Arabidopsis and Botrytis proteins
#This is the number of proteins for each organism in the raw data
df %>%
	distinct(protein_ID) %>%                 # keep only unique IDs
	mutate(prefix = str_sub(protein_ID, 1, 2)) %>%
	count(prefix)


df %>%
	write.csv("data/timecourse/input/AtBc_Proteome_TimeCourse_filtered.csv", row.names = F)
