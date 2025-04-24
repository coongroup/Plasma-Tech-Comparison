
library(devtools)
library(tidyverse)
library(RColorBrewer)



##* adjust column names 
seer_PG_NPs <- seer_PG_NPs %>% 
  add_column(R.Condition = 1, .after = 1)

canc1 <- list(neat = neat_PG, 
              acid = acid_PG, 
              magnet = magnet_PG, 
              preomics = preomics_PG, 
              seer = seer_PG_NPs)

for (i in names(canc1)) {
  
  colnames(canc1[[i]]) <- c("PG.ProteinGroups", "R.Condition", sapply(strsplit(colnames(canc1[[i]])[-c(1,2)], "_"), function(x) paste(x[5], collapse = "_")))
  
}

neat_PG <- canc1[["neat"]]
acid_PG <- canc1[["acid"]]
magnet_PG <- canc1[["magnet"]]
preomics_PG <- canc1[["preomics"]]
seer_PG_NPs <- canc1[["seer"]]


##* log2 transform 

canc1 <- list(neat = neat_PG, 
              acid = acid_PG, 
              magnet = magnet_PG, 
              preomics = preomics_PG, 
              seer = seer_PG_NPs)

for (i in names(canc1)) {
  
  canc1[[i]][,-c(1:2)] <- log2(canc1[[i]][,-c(1:2)])
  
}

neat_PG <- canc1[["neat"]]
acid_PG <- canc1[["acid"]]
magnet_PG <- canc1[["magnet"]]
preomics_PG <- canc1[["preomics"]]
seer_PG_NPs <- canc1[["seer"]]


##* Impute (ImpSeq) 

seer_PG_NPs <- seer_PG_NPs %>%
  mutate(R.Condition = as.integer(R.Condition))

canc1 <- list(neat = neat_PG, 
              acid = acid_PG, 
              magnet = magnet_PG, 
              preomics = preomics_PG, 
              seer = seer_PG_NPs,
              olink = olink_5_w)

for (i in names(canc1)) {
  
  w <- canc1[[i]]
  
  # Show how many non-NA values there are for each protein group in each study group
  w_na <- w %>%
    mutate(
      prop_NA_1_20 = rowMeans(!is.na(dplyr::select(., `1`:`20`))),
      prop_NA_21_40 = rowMeans(!is.na(dplyr::select(., `21`:`40`)))
    ) %>%
    mutate(max_NA_proportion = pmax(prop_NA_1_20, prop_NA_21_40)) %>%
    dplyr::select(-prop_NA_1_20, -prop_NA_21_40)
  
  # Make a list of IDs to keep where there are at least 50% non-NA values in one of the cohorts
  ids_to_keep <- w_na %>%
    filter(max_NA_proportion >= 0.5) %>% 
    pull(PG.ProteinGroups)
  
  w <- w %>%
    filter(PG.ProteinGroups %in% ids_to_keep)
  
  # impute with ImpSeq
  w1 <- w %>%
    dplyr::select(where(is.double))
  
  data_imp <- impSeq(w1)
  
  data_imp <- w %>%
    dplyr::select(PG.ProteinGroups) %>%
    cbind(data_imp)
  
  df_name <- paste0(i, "_PG_imp")
  assign(df_name, data_imp)
  
}


##* Olink Data Median Normalization 
olink_PG_imp <- olink_PG_imp %>%
  mutate(across(where(is.numeric), ~ . - median(., na.rm = TRUE)))




#### Figure S8-S13 ####
library(pheatmap)
library(gridExtra)
#get sample metadata correct
healthy_meta <- fread("P:/Projects/WFB_2024_SeerPlasmaTechComparison/CancervsHealthySamples/selected_cancer_matched_control_samples.csv",
                      data.table = F) %>%
  dplyr::select(`Lot Number`, `Inventory Barcode`, Age, Sex) %>%
  mutate(`Treatment Status` = "NA") %>%
  mutate(group = "healthy")
cancer_meta <- fread("P:/Projects/WFB_2024_SeerPlasmaTechComparison/CancervsHealthySamples/selected_cancer_samples.csv",
                     data.table = F) %>%
  dplyr::select(`Lot Number`, `Inventory Barcode`, Age, Sex, `Treatment Status`) %>%
  mutate(group = "cancer")
all_meta <- rbind(healthy_meta, cancer_meta)
manifest <- fread("P:/Projects/WFB_2024_SeerPlasmaTechComparison/CancervsHealthySamples/BarcodeManifest.csv",
                  data.table = F)
all_meta <- all_meta %>%
  left_join(manifest %>% dplyr::select(`Inventory Barcode`, `William Label`), by = "Inventory Barcode")
sample_annot <- all_meta %>%
  column_to_rownames(var = 'William Label') %>%
  dplyr::select(Age, Sex, `Treatment Status`, group) %>%
  mutate(Status = `Treatment Status`) %>%
  dplyr::select(-`Treatment Status`) %>%
  rename(Group = group) %>%
  mutate(Group = case_when(
    Group == "cancer" ~ "Cancer",
    Group == "healthy" ~ "Healthy",
    T ~ Group
  ))


# using imputed log2 tranformed normalized dataframes
canc <- list(neat = neat_PG_imp, 
             acid = acid_PG_imp, 
             preomics = preomics_PG_imp, 
             magnet = magnet_PG_imp, 
             seer = seer_PG_imp,
             olink = olink_PG_imp)

pheatmaps <- list()
for (i in names(canc)) {
  
  # choose a data frame
  ne <- canc[[i]] 
  
  pheats <- 
    pheatmap(ne[,2:41],
             #color = viridis(24, direction = 1, option = "plasma"),
             color = rev(colorRampPalette(brewer.pal(24, "RdYlBu"))(24)),
             cluster_rows = T,
             cluster_cols = T,
             treeheight_row = 0,
             treeheight_col = 10,
             show_rownames = F,
             show_colnames = T,
             scale = "row",
             annotation_col = sample_annot,
             fontsize = 12,
             fontsize_col = 12,
             border_color = NA,
             main = paste0(i),
             width = 4,
             height = 3)
  
  pheatmaps[[paste0("heatmap_unimputed_", i)]] <- pheats$gtable
  
}


for (i in 1:length(pheatmaps)) {
  
  heat <- pheatmaps[[i]]
  
  ggsave(paste0("reports/figures/Experiment5_All_Imputed_Normalized_Heatmap_", i, ".tiff"), 
         plot = heat, 
         width = 7, 
         height = 8,
         dpi = 300)
  
}


png("reports/figures/Experiment5_All_Imputed_Normalized_Heatmaps.png", width = 2000, height = 3000,
    units = "px")
grid.arrange(
  grobs = pheatmaps,
  nrow = 3)
dev.off()

