##### Botrytis Mutant Course Model
##### For modeling Botrytis protein abundance across multiple host genotypes
##### March 2026 AJM

### Define paths from arguments =====================
args <- commandArgs(trailingOnly = TRUE)

#assign input files to specific variables
bc_counts_file <- args[1]		# Counts file
output_dir <- args[2]     # Output directory

# Ensure output directory ends with a slash (for safe concatenation)
if (!grepl("/$", output_dir)) {
	output_dir <- paste0(output_dir, "/")
}

### Load packages and data ===================
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)

#load count file
message("Loading bc counts file: ", bc_counts_file)
bc_long <- read.csv(bc_counts_file)

### Reformat count data ======================

bc_long <- bc_long %>% 
	dplyr::rename(protein_ID = gene_ID) %>%
	dplyr::select(protein_ID, treatment, genotype, rep, abundance)

#filter to only infected samples
bc_long <- bc_long %>%
	filter(treatment == "infected")

bc <- bc_long %>%
	pivot_wider(names_from = protein_ID,
							values_from = abundance)

bc$treatment <- as.factor(bc$treatment)
bc$genotype <- as.factor(bc$genotype)
bc$rep <- as.factor(bc$rep)

# Convert columns that start with "B" to integers
bc <- bc %>%
	mutate(across(starts_with("B"), ~ as.numeric(.)))

### Define model function ======================

analyze_protein <- function(protein_ID, bc) {
	
	message("modeling ", protein_ID)
	
	# model
	formula <- as.formula(paste(
		protein_ID,
		"~ genotype"
	))
	
	model <- glmmTMB(formula, data = bc_filt, family = nbinom2) %>%
		suppressMessages()
	
	#collect convergence note
	conv_note <- model$fit$message
	
	## ---- ANOVA + variance ----
	anova <- car::Anova(model) %>%
		as.data.frame() %>%
		rownames_to_column("variable") %>%
		mutate(protein_ID = protein_ID) %>%
		mutate(convergence_note = conv_note)
	
	fixed_var <- diag(vcov(model)$cond) %>%
		as.data.frame() %>%
		rownames_to_column("term") %>%
		setNames(c("term", "variance"))
	
	var_sums <- fixed_var %>%
		summarise(
			genotype = sum(variance[grepl("^genotype", term) & !grepl(":", term)], na.rm = TRUE),
			intercept = sum(variance[grepl("Intercept", term)], na.rm = TRUE),
			tot_var = sum(genotype, intercept)
		) %>%
		mutate(across(-tot_var, ~ .x / tot_var)) %>%
		dplyr::select(-tot_var) %>%
		pivot_longer(everything(), names_to = "variable", values_to = "variance")
	
	anova <- full_join(anova, var_sums, by = "variable") %>%
		dplyr::select(protein_ID, everything())
	anova$protein_ID <- protein_ID
	
	## ---- EMMs ----
	## EMMs on original (log link) scale
	emm_log <- emmeans(model, ~ genotype) %>%
		summary() %>%
		as.data.frame() %>%
		dplyr::rename(emmean_log = emmean) %>%
		mutate(across(c(emmean_log, SE, asymp.LCL, asymp.UCL), ~ .x * 1e4),   # rescale
					 protein_ID = protein_ID)
	
	## EMMs on response (original) scale
	emm_resp <- emmeans(model, ~ genotype, type = "response") %>%
		summary() %>%
		as.data.frame() %>%
		dplyr::rename(emmean_response = response) %>%
		mutate(across(c(emmean_response, SE, asymp.LCL, asymp.UCL), ~ .x * 1e4),  # rescale
					 protein_ID = protein_ID)
	
	## ---- DEG ----
	#derived from link scale because we are getting fold changes.
	DEG <- emmeans(model, specs = "genotype") %>%
		contrast(method = "revpairwise") %>%
		summary() %>%
		mutate(
			log2FC = estimate / log(2),
			SE = SE / log(2),
			protein_ID = protein_ID
		) %>%
		dplyr::select(protein_ID, everything(), -estimate)
	
	list(
		anova = anova,
		emm_log = emm_log,
		emm_resp = emm_resp,
		DEG = DEG
	)
}


### Run function for Col0 vs tgg12 ===================================

#filter to only Col0 and tgg12
bc_filt <- bc %>%
	filter(genotype %in% c("Col0", "tgg12"))

# Remove any columns from bc_filt where all values are zero
bc_filt <- bc_filt[, !sapply(bc_filt, function(col) all(col == 0, na.rm = TRUE))]

#scale columns that start with "B", divide them each by 1e6.
bc_filt <- bc_filt %>%
	mutate(across(starts_with("B"), ~ .x / 1e4))

#get list of genes to run
proteins <- unique(bc_long$protein_ID)

proteins <- proteins[1:10] #subset for testing

#setup outputs
results <- list()
failed_proteins <- character()

#run
for (protein in proteins) {
	tryCatch(
			results[[protein]] <- analyze_protein(protein, bc_filt),
		error = function(e) {
			message("Error for protein ", protein, " — skipping")
			failed_proteins <<- c(failed_proteins, protein)
		}
	)
}

# Combine outputs
anova_all <- bind_rows(lapply(results, `[[`, "anova"))
emm_log_all   <- bind_rows(lapply(results, `[[`, "emm_log"))
emm_resp_all   <- bind_rows(lapply(results, `[[`, "emm_resp"))
DEG_all   <- bind_rows(lapply(results, `[[`, "DEG"))

#Do FDR correction (BH)
# Split data by variable type
anova_split <- split(anova_all, anova_all$variable)
# Apply FDR correction to each variable group
anova_fdr <- lapply(anova_split, function(x) {
  p_values <- x$`Pr(>Chisq)`
  x$p_adj <- p.adjust(p_values, method = "BH")
  return(x)})
# Recombine into single dataframe
anova_corrected <- do.call(rbind, anova_fdr)
# Reset row names
rownames(anova_corrected) <- NULL

#Replace DEG p value with anova genotype p value
# #remove DEG p value
# DEG_all <- subset(DEG_all, select = -p.value)
#get anova p values
anova_ps <- anova_corrected %>% filter(variable == "genotype") %>%
	dplyr::select(protein_ID, p_adj)
#join with bc anova
DEG_all <- left_join(DEG_all, anova_ps, by = "protein_ID")

#Put failed genes in a df
failed_proteins_df <- data.frame(failed_protein = failed_proteins, stringsAsFactors = FALSE)

#write out results
message("Writing results to output directory: ", paste0(output_dir, "bc_tgg12vsCol0/"))
dir.create(output_dir, recursive = T)
dir.create(paste0(output_dir, "bc_tgg12vsCol0/"), recursive = T)
write.csv(emm_log_all, paste0(output_dir, "bc_tgg12vsCol0/", "bc_tgg12vsCol0_adjusted_emmeans_log.csv"), row.names = F)
write.csv(emm_resp_all, paste0(output_dir, "bc_tgg12vsCol0/", "bc_tgg12vsCol0_adjusted_emmeans_response.csv"), row.names = F)
write.csv(anova_corrected, paste0(output_dir, "bc_tgg12vsCol0/", "bc_tgg12vsCol0_anova.csv"), row.names = F)
write.csv(DEG_all, paste0(output_dir, "bc_tgg12vsCol0/", "bc_tgg12vsCol0_DEGs.csv"), row.names = F)
write.csv(failed_proteins_df, paste0(output_dir, "bc_tgg12vsCol0/", "failed_proteins.csv"), row.names = FALSE)


### Run function for Col0 vs aop2 ===================================

#filter to only Col0 and aop2
bc_filt <- bc %>%
	filter(genotype %in% c("Col0", "aop2"))

# Remove any columns from bc_filt where all values are zero
bc_filt <- bc_filt[, !sapply(bc_filt, function(col) all(col == 0, na.rm = TRUE))]

#scale columns that start with "B", divide them each by 1e6.
bc_filt <- bc_filt %>%
	mutate(across(starts_with("B"), ~ .x / 1e4))

#get list of genes to run
proteins <- unique(bc_long$protein_ID)

proteins <- proteins[1:10] #subset for testing

#setup outputs
results <- list()
failed_proteins <- character()

#run
for (protein in proteins) {
	tryCatch(
		results[[protein]] <- analyze_protein(protein, bc_filt),
		error = function(e) {
			message("Error for protein ", protein, " — skipping")
			failed_proteins <<- c(failed_proteins, protein)
		}
	)
}

# Combine outputs
anova_all <- bind_rows(lapply(results, `[[`, "anova"))
emm_log_all   <- bind_rows(lapply(results, `[[`, "emm_log"))
emm_resp_all   <- bind_rows(lapply(results, `[[`, "emm_resp"))
DEG_all   <- bind_rows(lapply(results, `[[`, "DEG"))

#Do FDR correction (BH)
# Split data by variable type
anova_split <- split(anova_all, anova_all$variable)
# Apply FDR correction to each variable group
anova_fdr <- lapply(anova_split, function(x) {
	p_values <- x$`Pr(>Chisq)`
	x$p_adj <- p.adjust(p_values, method = "BH")
	return(x)})
# Recombine into single dataframe
anova_corrected <- do.call(rbind, anova_fdr)
# Reset row names
rownames(anova_corrected) <- NULL

#Replace DEG p value with anova genotype p value
# #remove DEG p value
# DEG_all <- subset(DEG_all, select = -p.value)
#get anova p values
anova_ps <- anova_corrected %>% filter(variable == "genotype") %>%
	dplyr::select(protein_ID, p_adj)
#join with bc anova
DEG_all <- left_join(DEG_all, anova_ps, by = "protein_ID")

#Put failed genes in a df
failed_proteins_df <- data.frame(failed_protein = failed_proteins, stringsAsFactors = FALSE)

#write out results
message("Writing results to output directory: ", paste0(output_dir, "bc_aop2vsCol0/"))
dir.create(paste0(output_dir, "bc_aop2vsCol0/"), recursive = T)
write.csv(emm_log_all, paste0(output_dir, "bc_aop2vsCol0/", "bc_aop2vsCol0_adjusted_emmeans_log.csv"), row.names = F)
write.csv(emm_resp_all, paste0(output_dir, "bc_aop2vsCol0/", "bc_aop2vsCol0_adjusted_emmeans_response.csv"), row.names = F)
write.csv(anova_corrected, paste0(output_dir, "bc_aop2vsCol0/", "bc_aop2vsCol0_anova.csv"), row.names = F)
write.csv(DEG_all, paste0(output_dir, "bc_aop2vsCol0/", "bc_aop2vsCol0_DEGs.csv"), row.names = F)
write.csv(failed_proteins_df, paste0(output_dir, "bc_aop2vsCol0/", "failed_proteins.csv"), row.names = FALSE)

### Run function for Col0 vs mybs ===================================

#filter to only Col0 and mybs
bc_filt <- bc %>%
	filter(genotype %in% c("Col0", "myb28", "myb29", "myb2829"))

# Remove any columns from bc_filt where all values are zero
bc_filt <- bc_filt[, !sapply(bc_filt, function(col) all(col == 0, na.rm = TRUE))]

#scale columns that start with "B", divide them each by 1e6.
bc_filt <- bc_filt %>%
	mutate(across(starts_with("B"), ~ .x / 1e4))

#get list of genes to run
proteins <- unique(bc_long$protein_ID)

proteins <- proteins[1:10] #subset for testing

#setup outputs
results <- list()
failed_proteins <- character()

#run
for (protein in proteins) {
	tryCatch(
		results[[protein]] <- analyze_protein(protein, bc_filt),
		error = function(e) {
			message("Error for protein ", protein, " — skipping")
			failed_proteins <<- c(failed_proteins, protein)
		}
	)
}

# Combine outputs
anova_all <- bind_rows(lapply(results, `[[`, "anova"))
emm_log_all   <- bind_rows(lapply(results, `[[`, "emm_log"))
emm_resp_all   <- bind_rows(lapply(results, `[[`, "emm_resp"))
DEG_all   <- bind_rows(lapply(results, `[[`, "DEG"))

#Do FDR correction (BH)
# Split data by variable type
anova_split <- split(anova_all, anova_all$variable)
# Apply FDR correction to each variable group
anova_fdr <- lapply(anova_split, function(x) {
	p_values <- x$`Pr(>Chisq)`
	x$p_adj <- p.adjust(p_values, method = "BH")
	return(x)})
# Recombine into single dataframe
anova_corrected <- do.call(rbind, anova_fdr)
# Reset row names
rownames(anova_corrected) <- NULL

#Replace DEG p value with anova genotype p value
# #remove DEG p value
# DEG_all <- subset(DEG_all, select = -p.value)
#get anova p values
anova_ps <- anova_corrected %>% filter(variable == "genotype") %>%
	dplyr::select(protein_ID, p_adj)
#join with bc anova
DEG_all <- left_join(DEG_all, anova_ps, by = "protein_ID")

#Put failed genes in a df
failed_proteins_df <- data.frame(failed_protein = failed_proteins, stringsAsFactors = FALSE)

#write out results
message("Writing results to output directory: ", paste0(output_dir, "bc_mybsvsCol0/"))
dir.create(paste0(output_dir, "bc_mybsvsCol0/"), recursive = T)
write.csv(emm_log_all, paste0(output_dir, "bc_mybsvsCol0/", "bc_mybsvsCol0_adjusted_emmeans_log.csv"), row.names = F)
write.csv(emm_resp_all, paste0(output_dir, "bc_mybsvsCol0/", "bc_mybsvsCol0_adjusted_emmeans_response.csv"), row.names = F)
write.csv(anova_corrected, paste0(output_dir, "bc_mybsvsCol0/", "bc_mybsvsCol0_anova.csv"), row.names = F)
write.csv(DEG_all, paste0(output_dir, "bc_mybsvsCol0/", "bc_mybsvsCol0_DEGs.csv"), row.names = F)
write.csv(failed_proteins_df, paste0(output_dir, "bc_mybsvsCol0/", "failed_proteins.csv"), row.names = FALSE)

### Run function for WT accessions ===================================

#filter to only WT accessions
bc_filt <- bc %>%
	filter(genotype %in% c("Col0","Ler","Sha"))

# Remove any columns from bc_filt where all values are zero
bc_filt <- bc_filt[, !sapply(bc_filt, function(col) all(col == 0, na.rm = TRUE))]

#scale columns that start with "B", divide them each by 1e6.
bc_filt <- bc_filt %>%
	mutate(across(starts_with("B"), ~ .x / 1e4))

#get list of genes to run
proteins <- unique(bc_long$protein_ID)

proteins <- proteins[1:10] #subset for testing

#setup outputs
results <- list()
failed_proteins <- character()

#run
for (protein in proteins) {
	tryCatch(
		results[[protein]] <- analyze_protein(protein, bc_filt),
		error = function(e) {
			message("Error for protein ", protein, " — skipping")
			failed_proteins <<- c(failed_proteins, protein)
		}
	)
}

# Combine outputs
anova_all <- bind_rows(lapply(results, `[[`, "anova"))
emm_log_all   <- bind_rows(lapply(results, `[[`, "emm_log"))
emm_resp_all   <- bind_rows(lapply(results, `[[`, "emm_resp"))
DEG_all   <- bind_rows(lapply(results, `[[`, "DEG"))

#Do FDR correction (BH)
# Split data by variable type
anova_split <- split(anova_all, anova_all$variable)
# Apply FDR correction to each variable group
anova_fdr <- lapply(anova_split, function(x) {
	p_values <- x$`Pr(>Chisq)`
	x$p_adj <- p.adjust(p_values, method = "BH")
	return(x)})
# Recombine into single dataframe
anova_corrected <- do.call(rbind, anova_fdr)
# Reset row names
rownames(anova_corrected) <- NULL

#Replace DEG p value with anova genotype p value
# #remove DEG p value
# DEG_all <- subset(DEG_all, select = -p.value)
#get anova p values
anova_ps <- anova_corrected %>% filter(variable == "genotype") %>%
	dplyr::select(protein_ID, p_adj)
#join with bc anova
DEG_all <- left_join(DEG_all, anova_ps, by = "protein_ID")

#Put failed genes in a df
failed_proteins_df <- data.frame(failed_protein = failed_proteins, stringsAsFactors = FALSE)

#write out results
message("Writing results to output directory: ", paste0(output_dir, "bc_accessions/"))
dir.create(paste0(output_dir, "bc_accessions/"), recursive = T)
write.csv(emm_log_all, paste0(output_dir, "bc_accessions/", "bc_accessions_adjusted_emmeans_log.csv"), row.names = F)
write.csv(emm_resp_all, paste0(output_dir, "bc_accessions/", "bc_accessions_adjusted_emmeans_response.csv"), row.names = F)
write.csv(anova_corrected, paste0(output_dir, "bc_accessions/", "bc_accessions_anova.csv"), row.names = F)
write.csv(DEG_all, paste0(output_dir, "bc_accessions/", "bc_accessions_DEGs.csv"), row.names = F)
write.csv(failed_proteins_df, paste0(output_dir, "bc_accessions/", "failed_proteins.csv"), row.names = FALSE)