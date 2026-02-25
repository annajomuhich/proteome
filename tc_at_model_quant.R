##### Arabidopsis Time Course Model
##### For modeling arabidopsis protein abundance over time
##### February 2026 AJM

### Define paths from arguments =====================
args <- commandArgs(trailingOnly = TRUE)

#assign input files to specific variables
at_counts_file <- args[1]		# Counts file (at_norm_counts_expressed.csv)
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

#load count file for arabidopsis
message("Loading at counts file: ", at_counts_file)
at_long <- read.csv(at_counts_file)

### Reformat count data ======================

at_long <- at_long %>%
	filter(organism == "arabidopsis",
	genotype == "Col.0")

at_long <- at_long %>% dplyr::select(protein_ID, treatment, hpi, rep, abundance)

at <- at_long %>%
	pivot_wider(names_from = protein_ID,
							values_from = abundance)

at$treatment <- as.factor(at$treatment)
at$hpi <- as.factor(at$hpi)
at$rep <- as.factor(at$rep)

# Convert columns that start with "AT" to integers
at <- at %>%
	mutate(across(starts_with("AT"), ~ as.integer(.)))


### Define model function ======================

analyze_protein <- function(protein_ID, at) {
	
	message("modeling ", protein_ID)
	
	# model
	formula <- as.formula(paste(
		protein_ID,
		"~ treatment + hpi + treatment * hpi"
	))
	
	model <- glmmTMB(formula, data = at, family = nbinom2) %>%
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
			treatment = sum(variance[grepl("^treatment", term) & !grepl(":", term)], na.rm = TRUE),
			hpi = sum(variance[grepl("^hpi", term)], na.rm = TRUE),
			`treatment:hpi` = sum(variance[grepl("^treatment", term) & grepl(":", term)], na.rm = TRUE),
			intercept = sum(variance[grepl("Intercept", term)], na.rm = TRUE),
			tot_var = sum(treatment, hpi, `treatment:hpi`, intercept)
		) %>%
		mutate(across(-tot_var, ~ .x / tot_var)) %>%
		dplyr::select(-tot_var) %>%
		pivot_longer(everything(), names_to = "variable", values_to = "variance")
	
	anova <- full_join(anova, var_sums, by = "variable") %>%
		dplyr::select(protein_ID, everything())
	anova$protein_ID <- protein_ID
	
	## ---- EMMs ----
	emm_log <- emmeans(model, ~ treatment * hpi) %>%
		summary() %>%
		as.data.frame() %>%
		dplyr::rename(emmean_log = emmean) %>%
		mutate(protein_ID = protein_ID)
	emm_resp <- emmeans(model, ~ treatment + hpi, type = "response") %>%
		summary() %>%
		as.data.frame() %>%
		dplyr::rename(emmean_response = response) %>%
		mutate(protein_ID = protein_ID)
	
	## ---- DEG ----
	#derived from link scale because we are getting fold changes.
	DEG <- emmeans(model, specs = "treatment") %>%
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


### Run function across genes ===================================

#get list of genes to run
proteins <- unique(at_long$protein_ID)

proteins <- proteins[1:10] #subset for testing

#setup outputs
results <- list()
failed_proteins <- character()

#run
for (protein in proteins) {
	tryCatch(
		results[[protein]] <- analyze_protein(protein, at),
		error = function(e) {
			message("Error for protein ", protein, " — skipping")
			failed_proteins <<- c(failed_proteins, protein)
		}
	)
}

### Reformat and write outputs ===========================

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
#remove DEG p value
DEG_all <- subset(DEG_all, select = -p.value)
#get anova p values
anova_ps <- anova_corrected %>% filter(variable == "treatment") %>%
	dplyr::select(protein_ID, p_adj)
#join with at anova
DEG_all <- left_join(DEG_all, anova_ps, by = "protein_ID")

#Put failed genes in a df
failed_proteins_df <- data.frame(failed_protein = failed_proteins, stringsAsFactors = FALSE)

#write out results
message("Writing results to output directory: ", output_dir)
dir.create(output_dir, recursive = T)
write.csv(emm_log_all, paste0(output_dir, "at_tc_adjusted_emmeans_log.csv"), row.names = F)
write.csv(emm_resp_all, paste0(output_dir, "at_tc_adjusted_emmeans_response.csv"), row.names = F)
write.csv(anova_corrected, paste0(output_dir, "at_tc_anova.csv"), row.names = F)
write.csv(DEG_all, paste0(output_dir, "at_tc_DEGs.csv"), row.names = F)
write.csv(failed_proteins_df, paste0(output_dir, "failed_proteins.csv"), row.names = FALSE)
