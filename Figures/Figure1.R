library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(OlinkAnalyze)
library(eulerr)



pal <- c("#23B56D",
         "#F26439",
         "#5762AC",
         "#E38CBB",
         "#2BA7DF",
         "#99C13C")



#### read in reports ####
neat_PG <- read.csv("data/processed/PG_Matrix_Experiment1_neat.csv")
acid_PG <- read.csv("data/processed/PG_Matrix_Experiment1_acid.csv")
preomics_PG <- read.csv("data/processed/PG_Matrix_Experiment1_preomics.csv")
magnet_PG <- read.csv("data/processed/PG_Matrix_Experiment1_magnet.csv")
seer_PG <- read.csv("data/processed/PG_Matrix_Experiment1_NPs_separate.csv")



#### Figure 1a ####
dop <- list(neat_PG = neat_PG, 
            acid_PG = acid_PG, 
            magnet_PG = magnet_PG, 
            preomics_PG = preomics_PG, 
            NPA_PG = NPA_PG,
            NPB_PG = NPB_PG,
            seer_PG_NPs = seer_PG_NPs)

metrics <- data.frame(
  method = character(),
  Column = character(),
  Non_NA_Count = numeric()
)

for (i in names(dop)) {
  
  x <- dop[[i]]
  
  # get complete cases
  complete_rows <- x %>%
    filter(complete.cases(.)) %>%  # Keep only rows with no NAs
    nrow()
  
  # get num ids per run (and total)
  x1 <- x %>%
    summarize(across(everything(), ~ sum(!is.na(.)))) %>% 
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Non_NA_Count") %>%
    bind_rows(tibble(Column = "complete", Non_NA_Count = complete_rows)) %>%
    mutate(method = i)
  
  metrics <- rbind(metrics, x1)
  
}


# add Seer combined and Olink
olink <- read_NPX("data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet")
olink <- olink_lod(olink, lod_method = "NCLOD")
olink_1 <- olink %>% 
  filter(grepl('O1', SampleID)) %>%
  mutate(olink_suggested_signal = PCNormalizedNPX > PCNormalizedLOD) %>%
  filter(!(AssayQC == "WARN" | AssayQC == "NA" | SampleQC == "FAIL")) %>%
  filter(olink_suggested_signal == T) %>%
  mutate(untransformedNPX = (2^PCNormalizedNPX))

olink_w <- olink_1 %>%
  dplyr::select(SampleID, UniProt, untransformedNPX) %>%
  pivot_wider(names_from = SampleID, values_from = untransformedNPX, values_fill = NA)
x <- olink_w 
complete_rows <- x %>%
  filter(complete.cases(.)) %>%  # Keep only rows with no NAs
  nrow()
x1 <- x %>%
  rename(PG.ProteinGroups = UniProt) %>%
  summarize(across(everything(), ~ sum(!is.na(.)))) %>% 
  pivot_longer(cols = everything(), names_to = "Column", values_to = "Non_NA_Count") %>%
  bind_rows(tibble(Column = "complete", Non_NA_Count = complete_rows)) %>%
  mutate(method = "olink_PG")
metrics <- rbind(metrics, x1)



tot_comp <- metrics %>%
  filter(Column %in% c("complete", "PG.ProteinGroups"))

metrics <- metrics %>%
  filter(!Column %in% c("complete", "PG.ProteinGroups"))


tot_comp$method <- factor(tot_comp$method, levels = c("neat_PG", "acid_PG", "preomics_PG", "magnet_PG", "NPA_PG", "NPB_PG", "seer_PG_NPs", "olink_PG"))
metrics$method <- factor(metrics$method, levels = c("neat_PG", "acid_PG", "preomics_PG", "magnet_PG", "NPA_PG", "NPB_PG", "seer_PG_NPs", "olink_PG"))

means <- metrics %>%
  group_by(method) %>%
  summarise(mean_value = mean(Non_NA_Count))

tot_comp1 <- tot_comp %>%
  pivot_wider(names_from = Column, values_from = Non_NA_Count)

ggplot(metrics %>% filter(method != "NPA_PG", method != "NPB_PG"), 
       aes(method, Non_NA_Count, fill = method)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           width = 0.8,
           alpha = 0.75) +  
  geom_point(alpha = 0.5,
             shape = 21,
             color = "black",
             stroke = 0.1,
             size = 1,
             position = position_jitterdodge(dodge.width = 1, jitter.width = 1)) +  
  geom_errorbar(data = tot_comp1 %>% filter(method != "NPA_PG", method != "NPB_PG"), 
                aes(x = method, y = complete, ymin = complete, ymax = complete), 
                width = 0.4, 
                color = "black",
                size = 0.1) +
  geom_errorbar(data = tot_comp1 %>% filter(method != "NPA_PG", method != "NPB_PG"), 
                aes(x = method, y = PG.ProteinGroups, ymin = PG.ProteinGroups, ymax = PG.ProteinGroups), 
                width = 0.4, 
                color = "black",
                size = 0.1) +
  geom_text(data = means %>% filter(method != "NPA_PG", method != "NPB_PG"),
            aes(method, mean_value, label = round(mean_value, 0)), 
            angle = 0, 
            color = "black",
            vjust = -1, 
            hjust = 0.5,
            size = 2) +
  labs(x = NULL,
       y = "Protein Groups") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")
  ) 
ggsave("reports/figures/Experiment1_ProteinGroupCounts_Average_SeerCombined_noNPANPB_olinklod_small.pdf", 
       width = 8, height = 4, units = "cm")



#### Figure 1b ####
m <- list(neat = neat_PG, 
          acid = acid_PG, 
          magnet = magnet_PG, 
          preomics = preomics_PG, 
          seer = seer_PG)

for (i in names(m)) {
  
  # remove duplicate rows
  spec_raw <- m[[i]][-1] %>%
    #select(-R.Condition) %>%
    pivot_longer(cols = -PG.ProteinGroups,
                 names_to = "file",
                 values_to = "abundance") %>%
    mutate(method = str_split(file, "_", simplify = TRUE)[, 4])
  
  # make separate dataframe
  df_name <- paste0(i, "_PG")
  assign(df_name, spec_raw)
  
  # save csv
  write.csv(spec_raw, file = paste0("data/processed/PG_Long_Experiment1_MethodComparison_", i, ".csv"))
  
}

## Read in Olink 
olink <- read_NPX("data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet")
olink <- olink_lod(olink, lod_method = "NCLOD")
olink_1 <- olink %>% 
  filter(grepl('O1', SampleID)) %>%
  mutate(olink_suggested_signal = PCNormalizedNPX > PCNormalizedLOD) %>%
  filter(!(AssayQC == "WARN" | AssayQC == "NA" | SampleQC == "FAIL")) %>%
  filter(olink_suggested_signal == T) %>%
  mutate(untransformedNPX = (2^PCNormalizedNPX))


# take each MS method and Olink, and expand the protein groups to not miss anything

f <- list(neat = neat_PG, 
          acid = acid_PG, 
          magnet = magnet_PG, 
          preomics = preomics_PG, 
          seer = seer_PG)

overlap <- list()

for (i in names(f)) {
  g <- f[[i]]
  h <- unique(g$PG.ProteinGroups)
  overlap[[i]] <- h
}

g <- olink_1
h <- unique(g$UniProt)
overlap <- c(overlap, list(olink = h))

## seer vs olink venn ##
set1 <- overlap$seer
set2 <- overlap$olink

list_of_vectors <- list(
  Seer = set1,
  Olink = set2
)

euler_plot <- euler(list_of_vectors)
pdf("reports/figures/Seer_Olink_Overlap_olinklod.pdf", width = 2, height = 2)
plot(euler_plot,
     quantities = T,
     fills = c("#8DA0CB", "#A6CEE3", col[8]))
dev.off()


## olink vs all ms methods ##
# make full MS method vector
sel_v <- overlap[names(overlap) != "olink"]
res <- unique(unlist(sel_v))

set1 <- res
set2 <- overlap$olink

list_of_vectors <- list(
  MS = set1,
  Olink = set2
)

euler_plot <- euler(list_of_vectors)
pdf("reports/figures/AllMSMethods_Olink_Overlap_olinklod.pdf", width = 2, height = 2)
plot(euler_plot,
     quantities = T,
     fills = c("#B2DF8A", "#A6CEE3", "#B3B3B3"))
dev.off()



#### Figure 1c-d ####
# make a ranking plot with all methods using quant and proteins from 2021 HPPP
HPPP <- read.csv("data/metadata/2021_HPPP_Protein_List.csv") 


# 2. Add HPPP quant to each of the protein groups in spec_PG (olink + MS) to rank them
colnames(HPPP)[1] <- "PG.ProteinGroups"
HPPP <- HPPP %>%
  mutate(rank = dense_rank(dplyr::desc(max_abundance)))

# make a dataframe from the overlap vectors
spec_PG_rank <- bind_rows(lapply(names(overlap), function(name) {
  data_frame(ID = names(overlap[[name]]), method = name, PG.ProteinGroups = overlap[[name]])
}), .id = "source")

# Join with the abundance data
HPPP_rank <- spec_PG_rank %>%
  select(method, PG.ProteinGroups) %>%
  separate_rows(PG.ProteinGroups, sep = ";") %>%
  separate_rows(PG.ProteinGroups, sep = "_") %>%
  group_by(method) %>%
  left_join(HPPP, by = "PG.ProteinGroups") %>%
  select(method, PG.ProteinGroups, max_abundance, rank)



# count how many values there are from HPPP in there
sum(is.na(HPPP_rank$max_abundance)) #2816 (number keeps changing)

# Calculate medians for each method
median_order <- HPPP_rank %>%
  group_by(method) %>%
  summarize(median_value = median(max_abundance, na.rm = TRUE)) %>%
  arrange(desc(median_value)) %>%
  pull(method)
HPPP_rank$method <- factor(HPPP_rank$method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink"))

# count points for each method
countse <- HPPP_rank %>%
  group_by(method) %>%
  summarise(sum = sum(!is.na(max_abundance)))

medis <- HPPP_rank %>%
  group_by(method) %>%
  summarize(median_value = median(max_abundance, na.rm = TRUE))



# individual plots combined

ptplot <- ggplot(HPPP_rank, aes(rank, max_abundance, color = method)) +
  geom_point(size = 0.5,
             alpha = 0.1) +
  scale_color_manual(values = pal) +
  labs(x = NULL, y = "HPPP abundance") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "none") +
  facet_grid(method ~.) +
  geom_text(data = countse, 
            aes(x = 4200, y = 4, label = sum), 
            inherit.aes = FALSE,
            size = 2) +
  geom_text(data = countse, 
            aes(x = 4000, y = 4, label = method), 
            inherit.aes = FALSE,
            size = 2)


bpplot <- ggplot(HPPP_rank, aes(method, max_abundance, fill = method)) +
  geom_violin(alpha = 0.5, 
              width = 1,
              size = 0.2,
              color = NA) +
  geom_boxplot(width = 0.08, 
               alpha = 0.5,
               outliers = F,
               size = 0.2,
               color = "black") +
  #scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "HPPP abundance") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.position = "none") +
  geom_text(data = medis,
            aes(y = median_value, 
                label = round(median_value, 2)),
            size = 2,
            hjust = -0.8)


ptplot
ggsave("reports/figures/Experiment1_Ranking.pdf", 
       width = 4, height = 4, units = "in")

bpplot
ggsave("reports/figures/Experiment1_Distribution.pdf", 
       width = 4, height = 2, units = "in")


#### Figure 1e ####

neat_PG <- read.csv("data/processed/PG_Matrix_Experiment3_neat.csv")
acid_PG <- read.csv("data/processed/PG_Matrix_Experiment3_acid.csv")
preomics_PG <- read.csv("data/processed/PG_Matrix_Experiment3_preomics.csv")
magnet_PG <- read.csv("data/processed/PG_Matrix_Experiment3_magnet.csv")
seer_PG <- read.csv("data/processed/PG_Matrix_Experiment3_seer.csv") #might need to combine this again after reading in and add columns



dop <- list(neat = neat_PG, 
            acid = acid_PG, 
            magnet = magnet_PG, 
            preomics = preomics_PG,
            seer = seer_PG_NPs)

metrics <- data.frame(
  method = character(),
  Column = character(),
  Non_NA_Count = numeric()
)

for (i in names(dop)) {
  
  x <- dop[[i]]
  
  # get complete cases
  complete_rows <- x %>%
    filter(complete.cases(.)) %>%  # Keep only rows with no NAs
    nrow()
  
  # get num ids per run (and total)
  x1 <- x %>%
    summarize(across(everything(), ~ sum(!is.na(.)))) %>% 
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Non_NA_Count") %>%
    bind_rows(tibble(Column = "complete", Non_NA_Count = complete_rows)) %>%
    mutate(method = i)
  
  metrics <- rbind(metrics, x1)
  
}


# add Seer combined and Olink
olink <- read_NPX("data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet")
olink <- olink_lod(olink, lod_method = "NCLOD")
olink_3 <- olink %>% 
  filter(grepl('O3', SampleID)) %>%
  mutate(olink_suggested_signal = PCNormalizedNPX > PCNormalizedLOD) %>%
  filter(!(AssayQC == "WARN" | AssayQC == "NA" | SampleQC == "FAIL")) %>%
  filter(olink_suggested_signal == T) %>%
  mutate(untransformedNPX = (2^PCNormalizedNPX))

olink_w <- olink_3 %>%
  dplyr::select(SampleID, UniProt, untransformedNPX) %>%
  pivot_wider(names_from = SampleID, values_from = untransformedNPX, values_fill = NA)
x <- olink_w 
complete_rows <- x %>%
  filter(complete.cases(.)) %>%  # Keep only rows with no NAs
  nrow()
x1 <- x %>%
  rename(PG.ProteinGroups = UniProt) %>%
  summarize(across(everything(), ~ sum(!is.na(.)))) %>% 
  pivot_longer(cols = everything(), names_to = "Column", values_to = "Non_NA_Count") %>%
  bind_rows(tibble(Column = "complete", Non_NA_Count = complete_rows)) %>%
  mutate(method = "olink_PG")
metrics <- rbind(metrics, x1)



metrics1 <- metrics %>%
  filter(!Column %in% c("complete", "PG.ProteinGroups", "X", "R.Condition")) %>%
  mutate(lip = word(Column, 5, sep = "_")) %>%
  mutate(lip = gsub("[a-zA-Z]", "", lip)) %>%
  mutate(lip = case_when(
    Column == "O3-1" ~ "1",
    Column == "O3-2" ~ "2",
    Column == "O3-3" ~ "3",
    T ~ lip
  )) %>%
  mutate(method = factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink_PG")))


means <- metrics1 %>%
  group_by(method, lip) %>%
  summarise(mean_value = mean(Non_NA_Count))


ggplot(metrics1, 
       aes(method, Non_NA_Count, fill = method, alpha = lip)) +
  geom_bar(stat = "summary", 
           fun = "mean", 
           position = position_dodge(width = 0.8),
           width = 0.7) +  
  geom_point(shape = 21,
             color = "black",
             stroke = 0.1,
             size = 1,
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +  
  geom_text(data = means,
            aes(method, mean_value, label = round(mean_value, 0)), 
            angle = 45, 
            color = "black",
            vjust = -1, 
            hjust = 0,
            size = 2,
            position = position_dodge(width = 0.8)) +
  labs(x = NULL,
       y = "Protein Groups") +
  scale_fill_manual(values = pal) +
  scale_alpha_manual(values = c(1, 0.75, 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,5000)) +
  theme_classic() +
  guides(fill = "none") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")
  ) 






