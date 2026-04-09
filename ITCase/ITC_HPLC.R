### Reformatting and plotting raw HPLC output
### April 2026 AJM

library(tidyverse)

df_og <- read.csv("data/HPLC/raw/OUTPT.CSV")

df <- df_og %>%
    dplyr::select(Directory, X17.2, X17.9) %>%
    dplyr::rename(breakdown_product = X17.2,
                  PEITC = X17.9) %>%
    mutate(Directory = 
        ifelse(Directory == "013-P1-H12-Enz_129.D",
                            "013-P1-H12-Enz129_4.D", Directory)) %>% #fix one directory entry
    mutate(
        enzyme_conc = as.numeric(str_extract(Directory, "(?<=Enz)\\d+")),
        rep = as.numeric(str_extract(Directory, "(?<=_)[0-9]+(?=\\.D$)"))
    ) %>% #get enzyme concentration and replicate number
    # Replace NA values in breakdown product and PEITC with 0
    mutate(breakdown_product = ifelse(is.na(breakdown_product), 0, breakdown_product),
           PEITC = ifelse(is.na(PEITC), 0, PEITC)) %>%
    #calculate percent conversion
    mutate(percent_conversion = (breakdown_product / (breakdown_product + PEITC)) * 100) %>%
    #remove MeOH blank
    filter(Directory != "001-P1-G10-MeOH_blank.D") %>%
    #replace NA values in enzyme_conc with 0
    mutate(enzyme_conc = ifelse(is.na(enzyme_conc), 0, enzyme_conc))

#write reformatted dataframe
df %>% write.csv("data/HPLC/reformatted/ITC_HPLC_reformatted.csv", row.names = F)

#plot the percent conversion vs enzyme concentration
df %>%
    group_by(enzyme_conc) %>%
    summarise(mean_percent_conversion = mean(percent_conversion),
              sd_percent_conversion = sd(percent_conversion)) %>%
    ggplot(aes(x = enzyme_conc, y = mean_percent_conversion)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean_percent_conversion - sd_percent_conversion,
                      ymax = mean_percent_conversion + sd_percent_conversion),
                  width = 0.2) +
    geom_smooth(method = "lm", se = FALSE, color = "grey") +
    labs(x = "BcSaxA Concentration (µg/mL)", y = "PEITC Conversion (%)") +
    theme_classic()
ggsave("figures/ITCase/ITC_HPLC_percent_conversion.png", height = 3, width = 3,
       dpi = 2500)
