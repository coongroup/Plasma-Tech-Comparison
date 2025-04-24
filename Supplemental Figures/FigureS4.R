
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(eulerr)
library(OlinkAnalyze)
library(data.table)



col <- brewer.pal(8, "Set2") 

# Manuscript Palette #
pal <- c("#23B56D",
         "#F26439",
         "#5762AC",
         "#E38CBB",
         "#2BA7DF",
         "#99C13C")



neat_PG <- read.csv("data/processed/PG_Matrix_Experiment1_neat.csv")
acid_PG <- read.csv("data/processed/PG_Matrix_Experiment1_acid.csv")
preomics_PG <- read.csv("data/processed/PG_Matrix_Experiment1_preomics.csv")
magnet_PG <- read.csv("data/processed/PG_Matrix_Experiment1_magnet.csv")
seer_PG <- read.csv("data/processed/PG_Matrix_Experiment1_NPs_separate.csv")

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


#### Figure S4 ####
HPPP <- read.csv("data/metadata/2021_HPPP_Protein_List.csv") 

# 1. Make a venn diagram of HPPP proteins vs the ones in my experiment (for both Olink vs HPPP and MS vs HPPP and all vs HPPP) ####
#Getting the vectors ready
set1 <- overlap$olink
set2 <- res 
set3 <- unique(unlist(overlap)) 
set4 <- HPPP$protein_identifier

list_of_vectors <- list(
  Olink = set1,
  MS = set2,
  All = set3,
  HPPP = set4
)

## Olink HPPP
euler_plot <- euler(list_of_vectors[c(1,4)])
pdf("reports/figures/HPPP_Olink_Overlap_olinkLOD.pdf")
plot(euler_plot,
     quantities = T,
     fills = c("#B2DF8A", "#A6CEE3", "#B3B3B3"))
dev.off()

## MS HPPP
euler_plot <- euler(list_of_vectors[c(2,4)])
pdf("reports/figures/HPPP_MS_Overlap_olinkLOD.pdf")
plot(euler_plot,
     quantities = T,
     fills = c("#E78AC3", "#A6CEE3", "#B3B3B3"))
dev.off()

## All HPPP
euler_plot <- euler(list_of_vectors[c(3,4)])
pdf("reports/figures/HPPP_All_Overlap_olinkLOD.pdf")
plot(euler_plot,
     quantities = T,
     fills = c("#FFD92F", "#A6CEE3", "#B3B3B3"))
dev.off()

