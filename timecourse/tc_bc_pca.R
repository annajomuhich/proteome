### Botrytis time course proteome PCA
### April 2026 - AJM

library(tidyverse)
library(gridExtra)
library(edgeR)
library(readxl)
library(RColorBrewer)
library(ggrepel)

# df <- read.csv("data/timecourse/tc_at_model_quant_20260225/at_tc_adjusted_emmeans_log.csv") %>%
# 	dplyr::select(treatment, hpi, emmean_log, protein_ID)

df <- read.csv("data/timecourse/input/AtBc_Proteome_TimeCourse_filtered.csv") %>%
	filter(genotype == "Col.0" & organism == "botrytis") %>%
	dplyr::select(protein_ID, treatment, hpi, rep, abundance)

#Check for any NAs (needs to be FALSE)
any(is.na(df))

#pivot wide
# df <- df %>%
# 	pivot_wider(names_from = protein_ID,
# 							values_from = emmean_log)
df <- df %>%
	pivot_wider(names_from = protein_ID,
							values_from = abundance)

#Check for any NAs (needs to be FALSE)
any(is.na(df))
colSums(is.na(df) > 0) #see by column


######### PCA

#get a numerical dataframe
# dfpca <- df[,-c(1,2)]
dfpca <- df[,-c(1,3)]

#run PCA
pca_result <- prcomp(dfpca, scale. = FALSE)

#take a look
pca_summary <- summary(pca_result)
pca_summary
#extract PC1 and PC2 proportions of variance for axis labels
pca_summary <- as.data.frame(pca_summary$importance)
pc1 <- pca_summary["Proportion of Variance", "PC1"]
pc2 <- pca_summary["Proportion of Variance", "PC2"]
pc1 <- pc1 * 100
pc2 <- pc2 * 100
pc1 <- signif(pc1, 3) #reduce to 3 significant figures
pc2 <- signif(pc2, 3)

#convert to dataframe
pca_df <- as.data.frame(pca_result$x)
#initial plot
ggplot(data = pca_df, aes(PC1, PC2)) + geom_point()

#append host geno and iso information
pca_df <- cbind(df[,c(1:3)], pca_df)

pca_df <- pca_df %>%
	mutate(hpi = as.factor(hpi),
				 treatment = as.factor(treatment))

#PC1 x PC2
pca_df %>% 
	ggplot(aes(PC1, PC2)) +
	geom_point(aes(color = hpi,
								 shape = treatment), size = 3) +
	xlab(paste0("PC1 - ", pc1, "%")) +
	ylab(paste0("PC2 - ", pc2, "%")) +
	theme_minimal() +
	scale_color_manual(values = c("red",
																"#f76f00",
																"#ebba34",
																"forest green",
																"turquoise",
																"blue",
																"#6d3a91")) +
	theme(
		legend.key.size = unit(0.4, "cm")
	) +
 labs(color = "Time (HPI)",
 		 shape = "Treatment")

ggsave("figures/timecourse/botrytis/tc_bc_pca.png", width = 4.5, height = 3)
