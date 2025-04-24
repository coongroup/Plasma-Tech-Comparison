
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(OlinkAnalyze)

col <- brewer.pal(8, "Set2")


# files
spec_PG <- read.csv("data/processed/PG_Matrix_Experiment1_separate_combined_noNPANPB.csv")
olink <- read_NPX("data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet")
contaminant <- read.csv("data/metadata/contamination_proteins.csv") %>%
  mutate(single_id = word(Protein_IDs, 1, sep = ";")) %>%
  mutate(single_id = word(single_id, 1, sep = "-"))

# REMOVE Replicate run FOR EACH METHOD 
spec_PG <- spec_PG %>% 
  pivot_longer(cols = -PG.ProteinGroups, names_to = "group", values_to = "value") %>%
  filter(grepl("_1|_2|_3|_4|_5", group)) %>%
  filter(!grepl("_0_2|_1_2", group)) %>%
  pivot_wider(names_from = group, values_from = value)


# check for exact matches to contamination vector
contaminant1 <- contaminant %>%
  mutate(is = single_id %in% spec_PG$PG.ProteinGroups) %>%
  mutate(single_id = case_when(
    single_id == "P0CG48" ~ "P0CG47;P0CG48;P62979;P62987",
    T ~ single_id
  )) %>% # this changes that id to the same as PG.ProteinGroups so there is a match
  mutate(is = single_id %in% spec_PG$PG.ProteinGroups) %>%
  filter(is == T) %>%
  mutate(PG.ProteinGroups = single_id) %>%
  select(PG.ProteinGroups, Type)


#### Figure S3 ####
# Organize long format data
contam_long <- spec_PG %>%
  #  left_join(olink_w, by = "PG.ProteinGroups") %>% # add olink
  pivot_longer(cols = -PG.ProteinGroups, names_to = "sample", values_to = "abundance") %>%
  mutate(method = word(sample, 4, sep = "_")) %>% # make method column
  #  mutate(method = case_when(
  #    grepl("O1-", sample ) ~ "olink",
  #    T ~ method
  #  )) %>%
  group_by(PG.ProteinGroups, method) %>% # get rid of NAs and mean of each set of 5 replicates
  summarise(mean_value = mean(log2(abundance), na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(mean_value)) %>%
  left_join(contaminant1, by = "PG.ProteinGroups") %>% # Add contaminant table
  mutate(Type = case_when(
    is.na(Type) ~ "No",
    T ~ Type
  )) %>%
  group_by(method) %>%
  mutate(rank = dense_rank(desc(mean_value))) %>% # ranks
  group_by(method) %>%
  arrange(mean_value) %>%
  mutate(percentile = rev(percent_rank(mean_value))) %>%
  mutate(Type = factor(Type, levels = c("Erythrocyte", "Platelet", "Coagulation", "Exosome", "Canonical", "No"))) %>%
  mutate(method = factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer")))

table(contam_long$Type) 




E <- ggplot(contam_long %>%
              filter(Type %in% c("No", "Erythrocyte")), 
            aes(percentile, mean_value, color = Type, size = Type, alpha = Type)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c(col[2], col[8])) +
  scale_size_manual(values = c(2, 0.5)) +
  scale_alpha_manual(values = c(0.5, 0.1)) +
  labs(x = NULL, 
       y = "abundance") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")
  ) +
  facet_grid(method ~ .) 

P <- ggplot(contam_long %>%
              filter(Type %in% c("No", "Platelet")), 
            aes(percentile, mean_value, color = Type, size = Type, alpha = Type)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c(col[3], col[8])) +
  scale_size_manual(values = c(2, 0.5)) +
  scale_alpha_manual(values = c(0.5, 0.1)) +
  labs(x = NULL, 
       y = NULL) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")
  ) +
  facet_grid(method ~ .) 

C <- ggplot(contam_long %>%
              filter(Type %in% c("No", "Coagulation")), 
            aes(percentile, mean_value, color = Type, size = Type, alpha = Type)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c(col[7], col[8])) +
  scale_size_manual(values = c(2, 0.5)) +
  scale_alpha_manual(values = c(0.5, 0.1)) +
  labs(x = NULL, 
       y = NULL) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")
  ) +
  facet_grid(method ~ .) 

Ex <- ggplot(contam_long %>%
               filter(Type %in% c("No", "Exosome")), 
             aes(percentile, mean_value, color = Type, size = Type, alpha = Type)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c(col[1], col[8])) +
  scale_size_manual(values = c(2, 0.5)) +
  scale_alpha_manual(values = c(0.5, 0.1)) +
  labs(x = NULL, 
       y = "abundance") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")
  ) +
  facet_grid(method ~ .) 

Ca <- ggplot(contam_long %>%
               filter(Type %in% c("No", "Canonical")), 
             aes(percentile, mean_value, color = Type, size = Type, alpha = Type)) +
  geom_point(shape = 16) +
  scale_color_manual(values = c(col[4], col[8])) +
  scale_size_manual(values = c(2, 0.5)) +
  scale_alpha_manual(values = c(0.5, 0.1)) +
  labs(x = NULL, 
       y = NULL) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 6),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")
  ) +
  facet_grid(method ~ .) 


E + P + C + Ex + Ca
