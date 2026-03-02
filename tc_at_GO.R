### GO enrichment analysis on host orthologs (Solanaceae)
### February 2026 AJM

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
library(data.table)
annot <- fread("data/gene_descriptions/ATH_GO_GOSLIM.txt", skip = 4)
annot <- annot %>%
	dplyr::select(V1, V6)
colnames(annot) <- c("gene", "GO")
#reformat GO from TAIR
df_annot <- annot %>%
	group_by(gene) %>%
	summarise(GO = paste(unique(GO), collapse = ","), .groups = "drop")

## setup geneID2GO for combined annotation
geneID2GO <- geneID2GO_setup(df_annot)

#get all genes
all_genes <- df_annot %>% pull(gene)


### Run GO enrichment on 120 high 48H cluster vs all genes  --------------------------------------------------------------------------

high_48hpi <- read.csv("data/genes_of_interest/high_48hpi_cluster.csv") %>% 
	mutate(protein_ID = sub("\\..*$", "", protein_ID)) %>%
	pull(protein_ID)

result <- run_topGO(high_48hpi, all_genes)
#filter to only significant results
result <- result %>% filter(p_adj < 0.05)
result %>% write.csv("data/timecourse/GO/at_high48hpi_GOEnriched.csv", row.names = F)



