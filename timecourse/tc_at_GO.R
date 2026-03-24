### GO enrichment analysis - Arabidopsis time course
### March 2026 AJM

library(tidyverse)
library(topGO)
library(GO.db)
library(data.table)

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
		
		results <- GenTable(GOdata,
												classicFisher = resultFisher,
												topNodes = topNodes)
		
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

annot <- fread("data/gene_descriptions/ATH_GO_GOSLIM.txt", skip = 4)
annot <- annot %>%
	dplyr::select(V1, V6)
colnames(annot) <- c("gene", "GO")
#reformat GO from TAIR
df_annot <- annot %>%
	group_by(gene) %>%
	summarise(GO = paste(unique(GO), collapse = ","), .groups = "drop")

#get all detected proteins in the dataset
detected_genes <- read.csv("data/timecourse/input/AtBc_Proteome_TimeCourse_filtered.csv") %>%
	filter(startsWith(protein_ID, "AT")) %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID)) %>%
	pull(protein_ID) %>%
	unique()
#filter annotation to detected proteins
df_annot <- df_annot %>%
	filter(gene %in% detected_genes)

## setup geneID2GO for combined annotation
geneID2GO <- geneID2GO_setup(df_annot)

### Run GO enrichment on all HCA clusters vs detected proteins ---------------------------------------------------------

#load cluster lists
cluster_df <- read.csv("data/timecourse/arabidopsis/tc_at_model_quant_20260311/at_tc_clusters.csv") %>%
	mutate(protein_ID = sub("\\..*$", "", protein_ID))

cluster_df <- cluster_df %>%
	mutate(cluster = as.character(cluster))

go_results_list <- list()

#Loop through each cluster and test GO enrichment
for (clust in unique(cluster_df$cluster)) {
	df_filtered <- cluster_df %>% filter(cluster == clust)
	cluster_genes <- df_filtered$protein_ID
	result <- run_topGO(cluster_genes, detected_genes)
	go_results_list[[as.character(clust)]] <- result
}

# Add 'cluster' as a column to each dataframe in go_results_list
go_results_annotated <- lapply(names(go_results_list), function(cluster) {
	df <- go_results_list[[cluster]]
	df$cluster <- cluster
	df
})
# Bind all annotated dataframes together
go_results_df <- dplyr::bind_rows(go_results_annotated)

#save all
go_results_df %>%
	write.csv("data/timecourse/arabidopsis/GO/tc_ac_GO_all.csv", row.names = F)

#save BP results
BPsig <- go_results_df %>%
	filter(Ontology == "BP" & p_adj < 0.05)

BPsig_list <- split(BPsig, BPsig$cluster)
# Save each dataframe in BPsig_list to a separate CSV file, named for its cluster
for (cluster in names(BPsig_list)) {
	write.csv(BPsig_list[[cluster]], file = paste0("data/timecourse/arabidopsis/GO/tc_ac_GO_BPsig_cluster", cluster, ".csv"), row.names = F)
}

#save MF results
MFsig <- go_results_df %>%
	filter(Ontology == "MF" & p_adj < 0.05)

MFsig_list <- split(MFsig, MFsig$cluster)
# Save each dataframe in BPsig_list to a separate CSV file, named for its cluster
for (cluster in names(BPsig_list)) {
	write.csv(MFsig_list[[cluster]], file = paste0("data/timecourse/arabidopsis/GO/MF/tc_ac_GO_MFsig_cluster", cluster, ".csv"), row.names = F)
}

### 20260301 - Run GO enrichment on 120 high 48H cluster vs all genes  --------------------------------------------------------------------------
# #get all genes
# all_genes <- df_annot %>% pull(gene)
# 
# high_48hpi <- read.csv("data/genes_of_interest/high_48hpi_cluster.csv") %>% 
# 	mutate(protein_ID = sub("\\..*$", "", protein_ID)) %>%
# 	pull(protein_ID)
# 
# result <- run_topGO(high_48hpi, all_genes)
# #filter to only significant results
# result <- result %>% filter(p_adj < 0.05)
# result %>% write.csv("data/timecourse/GO/at_high48hpi_GOEnriched.csv", row.names = F)



