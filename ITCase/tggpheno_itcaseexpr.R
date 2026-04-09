library(tidyverse)

pheno <- read.csv("data/Col0_tgg12_lesion_72HPI.csv") %>%
	dplyr::select(-c(Lesion_number, Lesion_size_pixels))

pheno <- pheno %>%
	group_by(Isolate_name, Host_genotype) %>%
	summarise(mean_lesion = mean(Lesion_size_mm2))

expr <- read.csv("../ITCase/ITCase_expression/Eudi_ITCase_expression.csv") %>%
	filter(species == "arabidopsis") %>%
	dplyr::select(-species)

df <- left_join(pheno, expr, join_by(Isolate_name == iso_name))

df <- df %>%
	filter(!is.na(ITCase_expression))

df$Host_genotype <- as.factor(df$Host_genotype)

model <- lm(mean_lesion ~ ITCase_expression * Host_genotype, data = df)
summary(model)

model <- lm(mean_lesion ~ Host_genotype, data = df)
summary(model)

df %>%
	ggplot(aes(x = ITCase_expression, y = mean_lesion, color = Host_genotype)) +
		geom_point() +
		geom_smooth(method = "lm", se = FALSE) +
		theme_minimal() +
		scale_color_manual(values = c("#81c97f", "#666666")) +
		labs(x = "ITCase Expression (log2)", y = expression("Mean Lesion Size (mm"^2*")"), color = "Host Genotype")

library(dplyr)
library(ggplot2)

# get R2 per genotype
r2_df <- df %>%
	group_by(Host_genotype) %>%
	summarise(
		model = list(lm(mean_lesion ~ ITCase_expression)),
		.groups = "drop"
	) %>%
	mutate(
		r2 = sapply(model, function(m) summary(m)$r.squared),
		p_value = sapply(model, function(m) summary(m)$coefficients[2, 4])
	) %>%
	select(-model)

# plot
library(ggrepel) # Make sure ggrepel is loaded for geom_text_repel

df %>%
	mutate(Host_genotype = recode(Host_genotype, "tgg12" = "tgg1/2")) %>%
	ggplot(aes(x = ITCase_expression, y = mean_lesion, color = Host_genotype)) +
	geom_point() +
	geom_smooth(method = "lm", se = FALSE) +
	# geom_text_repel(
	# 	data = . %>% filter(Isolate_name == "B05_10"),
	# 	aes(label = Isolate_name),
	# 	color = "black",
	# 	size = 3.5,
	# 	fontface = "bold",
	# 	inherit.aes = TRUE
	# ) +
	theme_classic() +
	scale_color_manual(values = c("Col0" = "#81c97f", "tgg1/2" = "#666666")) +
	labs(
		x = "ITCase Expression (log2)",
		y = expression("Mean Lesion Size (mm"^2*")"),
		color = "Host Genotype"
	)
ggsave("figures/ITCase/pheno_itcaseexpr.png", height = 4, width = 5)
