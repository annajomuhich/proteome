### Plotting Botrydial protein abundance vs ITCase protein abundance
### April 2026 AJM

df <- read.csv("data/mutant/input/Bc_Proteome_Mut_filtered.csv")

### ---------- Plot BcBot expression overall -----------------------------

df %>%
	mutate(genotype = factor(genotype, levels = c("Col0", "Ler", "Sha","AOP2", "myb28", "myb29", "myb2829", "tgg12"))) %>%
	filter(startsWith(Gene, "Bcbot")) %>%
	ggplot(aes(x = Gene, y = abundance, fill = genotype)) +
	geom_boxplot(size = 0.2) +   # reduced box outline thickness
	theme_classic() +
	labs(x = "", y = "Raw Protein Abundance") +
	scale_fill_brewer(palette = "Accent")
ggsave("figures/mutant/botrytis/mut_bc_botabundance.png", height = 4, width = 5,
       dpi = 2500)

### ---------- Plot BcBot expression vs ITCase expression -----------------------------

df_plot <- df %>%
	mutate(genotype = factor(genotype, levels = c("Col0", "Ler", "Sha","AOP2", "myb28", "myb29", "myb2829", "tgg12"))) %>%
	filter(startsWith(Gene, "Bcbot") | gene_ID == "Bcin06g00024") %>%
    mutate(Gene = if_else(gene_ID == "Bcin06g00024", "BcSaxA", Gene)) %>%
    dplyr::select(Gene, abundance, genotype, rep) %>%
    pivot_wider(names_from = Gene, values_from = abundance) %>%
    #scale each protein within each protein
    mutate(across(.cols = starts_with("Bc"), .fns = scale))

df_plot_long <- df_plot %>%
  pivot_longer(cols = starts_with("Bcbot"), 
               names_to = "Bcbot", 
               values_to = "Bcbot_abundance")

df_plot_long %>%
  ggplot(aes(x = Bcbot_abundance, y = BcSaxA, color = genotype)) +
  geom_point() +
  facet_wrap(
    ~ Bcbot, nrow = 2,
    labeller = as_labeller(function(x) gsub("^Bcbot:", "", x)),
    strip.position = "top"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),    # Removes box around facet titles
    strip.text = element_text(face = "bold"), # optionally bold text
    legend.position = c(0.85, 0.2), # position legend low and left
    legend.title = element_text(size = 8, face = "bold"), # reduce legend title size
    legend.text = element_text(size = 8),  # reduce legend label size
    legend.key.size = unit(0.8, "lines")  # shrink key size
  ) +
  labs(x = "Bcbot Abundance (scaled)", y = "BcSaxA Abundance (scaled)") +
  scale_color_brewer(palette = "Accent")

ggsave("figures/mutant/botrytis/mut_bc_botvsitcase.png",
			 width = 6, height = 4)
