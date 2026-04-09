### GO enrichment analysis on Botrytis genes
### January 2026 AJM

library(tidyverse)
library(topGO)
library(GO.db)

#define function to setup geneID2GO using the annotation file
geneID2GO_setup <- function(df_annot) {
	df_annot <- df_annot %>%
		dplyr::select(gene, GO) %>%
		distinct(gene, .keep_all = TRUE) %>%
		filter(GO != "")
	# Rename gene column
	#names(df_annot)[names(df_annot) == "locusName"] <- "gene"
	# Build geneID2GO list directly in memory
	geneID2GO <- setNames(
		strsplit(df_annot$GO, "\\s*,\\s*|\\s+"), #split string by comma or space
		df_annot$gene
	)
	return(geneID2GO)
}

#define function to run GO enrichment between two groups
#loops through each ontology and compiles results for each
run_topGO <- function(groupA, groupB, topNodes = 100) {
	
	# get entire list of genes to compare
	all_genes <- unique(c(groupA, groupB))
	
	# set up binary vector
	geneList <- factor(as.integer(all_genes %in% groupA))
	names(geneList) <- all_genes
	
	ontologies <- c("BP", "MF", "CC")
	
	results_list <- lapply(ontologies, function(ont) {

		message("Running GO enrichment for ", ont, "...")
		
		GOdata <- new("topGOdata",
									ontology = ont,
									allGenes = geneList,
									geneSel = function(p) p == 1,
									annot = annFUN.gene2GO,
									gene2GO = geneID2GO,
									nodeSize = 10)
		
		resultFisher <- runTest(GOdata,
														algorithm = "classic",
														statistic = "fisher")
		
		available_nodes <- length(score(resultFisher))
		
		results <- GenTable(GOdata,
												classicFisher = resultFisher,
												topNodes = min(topNodes, available_nodes))
		
		# full GO terms
		full_terms <- Term(GOTERM[results$GO.ID])
		results$Term <- full_terms
		
		# adjust p-values
		results$p_adj <- p.adjust(results$classicFisher, method = "BH")
		
		# add ontology column
		results$Ontology <- ont
		results <- results %>% dplyr::select(Ontology, everything())
		
		results
	})
	
	# combine all ontologies
	dplyr::bind_rows(results_list)
}

### Prepare gene universe (geneID2GO) -----------------------------------------------
#load annotation file
df_annot <- read.csv("data/gene_descriptions/Bcin_Annotations_Full_transcript.csv")
# Reformat annotation
df_annot <- df_annot %>%
	distinct(X.Gene.ID., .keep_all = TRUE) %>% # Removing duplicates based on 'gene'
	dplyr::select(X.Gene.ID.,X.Computed.GO.Component.IDs.,X.Computed.GO.Function.IDs.,X.Computed.GO.Process.IDs.) # Selecting only the 'gene' and 'GO' columns
# Rename the X.Gene.ID. column to gene
names(df_annot)[names(df_annot) == "X.Gene.ID."] <- "gene"
# Combine the two columns into one "GO" column, removing "N/A"
df_annot <- df_annot %>%
	mutate(
		GO = paste(
			ifelse(X.Computed.GO.Function.IDs. != "N/A", X.Computed.GO.Function.IDs., ""),
			ifelse(X.Computed.GO.Process.IDs. != "N/A", X.Computed.GO.Process.IDs., ""),
			ifelse(X.Computed.GO.Component.IDs. != "N/A", X.Computed.GO.Component.IDs., ""),
			sep = ",") %>%
			gsub("(^,|,$|,,)", "", .)) %>%  # Remove leading/trailing commas and duplicate commas
	dplyr::select(gene, GO)  # Keep only the Gene ID and the new GO column

#get all detected proteins in the dataset
detected_genes <- read.csv("data/timecourse/input/AtBc_Proteome_TimeCourse_filtered.csv") %>%
	filter(startsWith(protein_ID, "B")) %>%
	pull(protein_ID) %>%
	unique()
#filter annotation to detected proteins
df_annot <- df_annot %>%
	filter(gene %in% detected_genes)

#setup geneID2GO
geneID2GO <- geneID2GO_setup(df_annot)
#get all genes
all_genes <- df_annot %>% pull(gene)

### Run GO enrichment on each first detection time vs all Botrytis genes -------------------------------------------

df <- read.csv("data/timecourse/botrytis/tc_bc_model_timepoint_20260309/bc_first_hpi.csv")
go_results_list <- list()

for (hpi in unique(df$first_hpi)) {
	df_hpi <- df %>% filter(first_hpi == hpi)
	hpi_genes <- df_hpi$gene_ID
	result <- run_topGO(hpi_genes, all_genes)
	go_results_list[[as.character(hpi)]] <- result
}

# Add 'hpi' as a column to each dataframe in go_results_list
go_results_annotated <- lapply(names(go_results_list), function(hpi) {
	df <- go_results_list[[hpi]]
	df$hpi <- hpi
	df
})
# Bind all annotated dataframes together
go_results_df <- dplyr::bind_rows(go_results_annotated)

#save all
go_results_df %>%
	write.csv("data/timecourse/botrytis/GO/tc_bc_GO_all.csv", row.names = F)

BPsig <- go_results_df %>%
	filter(p_adj < 0.01 &
				 	Ontology == "BP")

BPsig_list <- split(BPsig, BPsig$hpi)
# Save each dataframe in BPsig_list to a separate CSV file, named for its hpi
for (hpi in names(BPsig_list)) {
	write.csv(BPsig_list[[hpi]], file = paste0("data/timecourse/botrytis/GO/tc_bc_GO_BPsig_hpi_", hpi, ".csv"), row.names = F)
}

#Trimming contents of each based on Revigo output from web app.

### ------------------- tested against all proteins ----------------------------------
#16hpi had no redundancies.
#24hpi:
revigo_24hpi <- read.delim("data/timecourse/botrytis/GO/Revigo_BP_Table_24hpi.tsv") %>% pull(TermID)
BPsig_list[["24"]] %>%
	filter(GO.ID %in% revigo_24hpi) %>%
	write.csv("data/timecourse/botrytis/GO/tc_bc_GO_BPsig_hpi_24_revigo.csv", row.names = F)
#32hpi:
revigo_32hpi <- read.delim("data/timecourse/botrytis/GO/Revigo_BP_Table_32hpi.tsv") %>% pull(TermID)
BPsig_list[["32"]] %>%
	filter(GO.ID %in% revigo_32hpi) %>%
	write.csv("data/timecourse/botrytis/GO/tc_bc_GO_BPsig_hpi_32_revigo.csv", row.names = F)
#40hpi:
revigo_40hpi <- read.delim("data/timecourse/botrytis/GO/Revigo_BP_Table_40hpi.tsv") %>% pull(TermID)
BPsig_list[["40"]] %>%
	filter(GO.ID %in% revigo_40hpi) %>%
	write.csv("data/timecourse/botrytis/GO/tc_bc_GO_BPsig_hpi_40_revigo.csv", row.names = F)
#8 and 48hpi had no significant GO terms

### ----------------- tested against only detected proteins -----------------------


