
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(OlinkAnalyze)
library(data.table)
library(viridis)
library(patchwork)

#Colors#
#Make a classic palette
col <- brewer.pal(8, "Set2") 



#### Figure S14a ####

olink <- read_NPX("data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet")
olink <- olink_lod(olink, lod_method = "NCLOD")
olink_5 <- olink %>% 
  filter(grepl('O5', SampleID)) %>%
  mutate(olink_suggested_signal = PCNormalizedNPX > PCNormalizedLOD) %>%
  filter(!(AssayQC == "WARN" | AssayQC == "NA" | SampleQC == "FAIL")) %>%
  filter(olink_suggested_signal == T) %>%
  mutate(untransformedNPX = (2^PCNormalizedNPX)) %>%
  mutate(
    extracted_num = as.numeric(str_extract(SampleID, "(?<=O5-)\\d+")),
    condition = case_when(
      extracted_num >= 1 & extracted_num <= 20 ~ "healthy",
      extracted_num >= 21 & extracted_num <= 40 ~ "cancer",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!(grepl('_r', SampleID)))


ol5 <- olink_5 %>%
  group_by(SampleID) %>%
  mutate(median_abundance = median(NPX, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(SampleID = factor(SampleID, levels = unique(SampleID[order(median_abundance)])))

p1 <- ggplot(ol5, aes(SampleID, NPX, fill = condition)) +
  geom_boxplot(width = 0.5, 
               alpha = 0.5, 
               outliers = F,
               size = 0.2) +
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        panel.spacing = unit(1, "lines"), 
        panel.grid.major.y = element_line(color = "gray", size = 0.5),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank()) 

# median adjustment
median_NPX <- olink_5 %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(Median_NPX = median(NPX)) 

# Adjust by sample median
olink_5_norm <- olink_5 %>% 
  dplyr::inner_join(median_NPX, by = "SampleID")|> 
  dplyr::mutate(NPX = NPX - Median_NPX) %>%
  group_by(SampleID) %>%
  mutate(median_abundance = median(NPX, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(SampleID = factor(SampleID, levels = unique(SampleID[order(median_abundance)])))

p2 <- ggplot(olink_5_norm, aes(SampleID, NPX, fill = condition)) +
  geom_boxplot(width = 0.5, 
               alpha = 0.5, 
               outliers = F,
               size = 0.2) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        panel.spacing = unit(1, "lines"), 
        panel.grid.major.y = element_line(color = "gray", size = 0.5),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank()) 

p1/p2