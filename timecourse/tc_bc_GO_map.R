### Presence absence matrix of Botrytis GO terms host specific and iso specific
### January 2026 AJM

library(tidyverse)

go_16 <- read.csv("data/timecourse/botrytis/GO/tc_bc_GO_BPsig_hpi_16_revigo.csv") %>%
	pull(Term)
go_24 <- read.csv("data/timecourse/botrytis/GO/tc_bc_GO_BPsig_hpi_24_revigo.csv") %>%
	pull(Term)
go_32 <- read.csv("data/timecourse/botrytis/GO/tc_bc_GO_BPsig_hpi_32_revigo.csv") %>%
	pull(Term)
go_40 <- read.csv("data/timecourse/botrytis/GO/tc_bc_GO_BPsig_hpi_40_revigo.csv") %>%
	pull(Term)
all_go <- sort(unique(c(go_16, go_24, go_32, go_40)))

#create empty vectors for 8 and 48
go_8 <- c()
go_48 <- c()
go_pa <- data.frame(GO = all_go,
										hpi8 = all_go %in% go_8,
										hpi16 = all_go %in% go_16,
										hpi24 = all_go %in% go_24,
										hpi32 = all_go %in% go_32,
										hpi40 = all_go %in% go_40,
										hpi48 = all_go %in% go_48)

go_long <- pivot_longer(
	go_pa,
	cols = -GO,
	names_to = "Timepoint",
	values_to = "Present"
)

go_long <- go_long %>%
	mutate(Present = if_else(Present == TRUE, "yes", "no")) %>%
	mutate(Timepoint = gsub("hpi", "", Timepoint))

# Set Timepoint as an ordered factor with correct numerical order
go_long <- go_long %>%
	mutate(Timepoint = factor(Timepoint, levels = c("8", "16", "24", "32", "40", "48"), ordered = TRUE))


ggplot(go_long, aes(x = Timepoint, y = GO, fill = Present)) +
	geom_tile(color = "grey80") +
	scale_fill_manual(
		values = c("no" = "white", "yes" = "dark orange"),
		name = "Upregulated"
	) +
	theme_minimal() +
	theme(
		axis.text.y = element_text(size = 7),
		panel.grid = element_blank(),
		legend.position = "none"
		#axis.text.x = element_text(angle = 45, hjust = 1)
	) +
	scale_y_discrete(limits = rev) +
	ylab("") +
	xlab("First Detection Time (hpi)")

ggsave("figures/timecourse/botrytis/tc_bc_GO_map.png", height = 6, width = 4.5)
