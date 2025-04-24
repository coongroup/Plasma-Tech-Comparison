
library(devtools)
library(tidyverse)
library(UpSetR)
library(OlinkAnalyze)
library(data.table)




#### Figure S2c ####

neat_PG <- fread("data/processed/IDmapping_Neat.tsv")
acid_PG <- fread("data/processed/IDmapping_Acid.tsv")
preomics_PG <- fread("data/processed/IDmapping_PreOmics.tsv")
magnet_PG <- fread("data/processed/IDmapping_magnet.tsv")
seer_PG <- fread("data/processed/IDmapping_Seer.tsv")
olink_PG <- fread("data/processed/IDmapping_Olink.tsv")

# take each MS method and Olink, and make a list

f <- list(neat = neat_PG, 
          acid = acid_PG, 
          magnet = magnet_PG, 
          preomics = preomics_PG, 
          seer = seer_PG,
          olink = olink_PG)

gene_overlap <- list()

for (i in names(f)) {
  g <- f[[i]] %>%
    dplyr::select(`Gene Names (primary)`)
  h <- unique(g$`Gene Names (primary)`)
  gene_overlap[[i]] <- h
}


pdf("reports/figures/Experiment1_average_method_PG_upset_small_genelevel.pdf", 
    width = 6.5, height = 3)
upset(fromList(gene_overlap), 
      sets = names(gene_overlap), 
      nintersects = NA,
      order.by = "freq",
      point.size = 1,
      line.size = 0.2,
      text.scale = 1,
      mainbar.y.label = 'Proteins',
      sets.x.label = 'Proteins per Method',
      number.angles = 45)

