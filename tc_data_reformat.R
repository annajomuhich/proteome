### Timecourse data reformat
### February 2026 AJM 

library(tidyverse)

df_raw <- read.delim("data/timecourse/0466_DIA-FA-GL_Botrytis-timecourse_54Samples_IntMatrix.tsv")

# Remove rows that are entirely NA except for column PG.ProteinGroups
df <- df_raw[!apply(df_raw[, setdiff(names(df_raw), "PG.ProteinGroups")], 1, function(x) all(is.na(x))), ]

# remove contaminant proteins
df <- df %>%
	filter(grepl("^AT", PG.ProteinGroups) |
				 	grepl("^XP", PG.ProteinGroups))

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
	select(protein_ID, PG.ProteinGroups, organism, sample_ID, genotype, treatment, hpi, rep, abundance)

#Replace NaN with 0
df <- df %>%
	mutate(abundance = if_else(is.na(abundance), 0, abundance))

### Low protein abundance filtering
#Get a list of proteins to keep. These are proteins that:
#For any combo of protein x treatment x hpi, there is signal in >= 2 reps
proteins_keep <- df %>%
	filter(abundance > 0) %>%
	group_by(protein_ID, treatment, hpi) %>%
	summarise(n_rep = n_distinct(rep), .groups = "drop") %>%
	filter(n_rep >= 2) %>%
	pull(protein_ID) %>%
	unique()

#remove low abundance proteins
df <- df%>%
	filter(protein_ID %in% proteins_keep)

df %>%
	write.csv("data/timecourse/AtBc_Proteome_TimeCourse_filtered.csv", row.names = F)
