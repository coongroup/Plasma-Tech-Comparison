
library(OlinkAnalyze)
require(tidyverse)
require(data.table)
library(ggrepel)
library(RColorBrewer)

pal <- c("#23B56D",
         "#F26439",
         "#5762AC",
         "#E38CBB",
         "#2BA7DF",
         "#99C13C")


methods <- c("neat", "acid", "preomics", "magnet", "NPA", "NPB")

# Make Correlation plots with concentration and protein quant for each method ----

MS_cor <- data.frame(Run = character(),
                     sample = character(),
                     Protein.Group = character(),
                     Protein.Names = character(),
                     log_abundance = numeric(),
                     nonLog2 = numeric(),
                     Group = character(),
                     method = character(),
                     concentration = numeric(),
                     cor = numeric())

for (i in methods) {
  
  rollupdata <- fread(paste0("data/diann_output/exp4/", i, "/report_filtered_withspecies_rollup.csv"))[,-1] %>%
    filter(str_detect(Protein.Names, "HUMAN")) %>%
    mutate(Group = substr(sample, 2,2)) %>%
    mutate(method = i,
           concentration = case_when(
             Group == 'a' ~ 100,
             Group == 'b' ~ 70,
             Group == 'c' ~ 50,
             Group == 'd' ~ 30,
             Group == 'e' ~ 10,
             Group == 'f' ~ 5,
             Group == 'g' ~ 1,
             Group == 'h' ~ .5,
             Group == 'i' ~ .1,
             Group == 'j' ~ 0
           )) 
  
  rollupdata <- rollupdata %>%
    group_by(Protein.Group) %>%
    mutate(cor = cor(concentration, nonLog2)) %>%
    ungroup
  
  MS_cor <- rbind(MS_cor, rollupdata)
  
}

# filter for NPA or NPB by more observations per protein group
resultasd <- MS_cor %>%
  filter(method %in% c("NPA", "NPB")) %>%
  group_by(Protein.Group) %>%
  filter(cor == max(cor)) 

MS_cor1 <- MS_cor %>%
  filter(!(method %in% c("NPA", "NPB"))) %>%
  bind_rows(resultasd) %>%
  mutate(method = case_when(
    method == "NPA" ~ "seer",
    method == "NPB" ~ "seer",
    TRUE ~ method
  ))


# Add Olink data (post-LOD calc)
dat <- read_parquet('data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet') 
dat <- olink_lod(dat, lod_method = "NCLOD")
dat_1 <- dat %>% 
  mutate(olink_suggested_signal = PCNormalizedNPX > PCNormalizedLOD) %>%
  filter(!(AssayQC == "WARN" | AssayQC == "NA" | SampleQC == "FAIL")) %>%
  filter(SampleType == 'SAMPLE') %>%
  filter(!grepl('_r', SampleID)) %>%
  filter(grepl('O4', SampleID)) %>%
  mutate(Group = substr(SampleID, 4,4)) %>%
  filter(olink_suggested_signal == T) %>%
  mutate(concentration = case_when(
    Group == 'a' ~ 100,
    Group == 'b' ~ 70,
    Group == 'c' ~ 50,
    Group == 'd' ~ 30,
    Group == 'e' ~ 10,
    Group == 'f' ~ 5,
    Group == 'g' ~ 1,
    Group == 'h' ~ .5,
    Group == 'i' ~ .1,
    Group == 'j' ~ 0
  )) %>%
  group_by(Assay) %>%
  mutate(cor = cor(concentration, PCNormalizedNPX)) %>%
  ungroup()

# make olink data to bind to MS_cor1
dat_1 <- dat_1 %>%
  dplyr::select(SampleID, UniProt, Assay, PCNormalizedNPX, Group, concentration, cor) %>%
  mutate(method = "olink",
         Run = SampleID,
         sample = SampleID,
         Protein.Group = UniProt,
         Protein.Names = Assay,
         log_abundance = PCNormalizedNPX,
         nonLog2 = (2^PCNormalizedNPX),
  ) %>%
  dplyr::select(-SampleID, -UniProt, -Assay, -PCNormalizedNPX)

MS_cor1 <- rbind(MS_cor1, dat_1)

MS_cor1$method <- factor(MS_cor1$method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink"))

#### Figure S7a ####
proportions <- MS_cor1 %>%
  group_by(method) %>%
  summarise(
    below_cutoff = mean(cor < 0.6, na.rm = T),
    above_cutoff = mean(cor >= 0.6, na.rm = T)
  )


ggplot(MS_cor1, aes(cor)) +
  geom_density() + 
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "red") +
  theme_bw() +
  facet_wrap(. ~ method) +
  geom_text(data = proportions,
            aes(x = 0.72, y = 4.5, label = round(above_cutoff, 2)),
            color = "black",
            size = 2) 



#### Figure S7b ####
dat <- read_parquet('data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet') 
dat <- olink_lod(dat, lod_method = "NCLOD")
dat_1 <- dat %>% 
  mutate(olink_suggested_signal = PCNormalizedNPX > PCNormalizedLOD) %>%
  filter(!(AssayQC == "WARN" | AssayQC == "NA" | SampleQC == "FAIL")) %>%
  filter(SampleType == 'SAMPLE') %>%
  filter(!grepl('_r', SampleID)) %>%
  filter(grepl('O4', SampleID)) %>%
  mutate(Group = substr(SampleID, 4,4)) %>%
  mutate(concentration = case_when(
    Group == 'a' ~ 100,
    Group == 'b' ~ 70,
    Group == 'c' ~ 50,
    Group == 'd' ~ 30,
    Group == 'e' ~ 10,
    Group == 'f' ~ 5,
    Group == 'g' ~ 1,
    Group == 'h' ~ .5,
    Group == 'i' ~ .1,
    Group == 'j' ~ 0
  )) %>%
  mutate(
    facet_group = case_when(
      olink_suggested_signal == TRUE  ~ "Above",
      olink_suggested_signal == FALSE ~ "Below",
      TRUE              ~ NA_character_
    )
  ) %>%
  bind_rows(mutate(., facet_group = "All")) %>%
  group_by(Assay) %>%
  mutate(cor = cor(concentration, PCNormalizedNPX)) %>%
  ungroup()


proportions1 <- dat_1 %>%
  group_by(facet_group) %>%
  summarise(
    below_cutoff = mean(cor < 0.6, na.rm = T),
    above_cutoff = mean(cor >= 0.6, na.rm = T)
  )



ggplot(dat_1, aes(cor)) +
  geom_density() + 
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "red") +
  theme_bw() +
  facet_wrap(. ~ facet_group) +
  geom_text(data = proportions1,
            aes(x = 0.72, y = 4.5, label = round(above_cutoff, 2)),
            color = "black",
            size = 2) 


