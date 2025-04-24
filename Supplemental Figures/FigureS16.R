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




neat_PG <- fread("data/processed/PG_Matrix_Experiment5_neat.csv", drop = c(1), data.table = F)
acid_PG <- fread("data/processed/PG_Matrix_Experiment5_acid.csv", drop = c(1), data.table = F)
preomics_PG <- fread("data/processed/PG_Matrix_Experiment5_preomics.csv", drop = c(1), data.table = F)
magnet_PG <- fread("data/processed/PG_Matrix_Experiment5_magnet.csv", drop = c(1), data.table = F)
seer_PG_NPs <- fread("data/processed/PG_Matrix_Experiment5_seer_PG_NPs.csv", drop = c(1), data.table = F)

olink <- read_NPX("data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet")
olink <- olink_lod(olink, lod_method = "NCLOD")
olink_5 <- olink %>% 
  filter(grepl('O5', SampleID)) %>%
  mutate(olink_suggested_signal = PCNormalizedNPX > PCNormalizedLOD) %>%
  filter(!(AssayQC == "WARN" | AssayQC == "NA" | SampleQC == "FAIL")) %>%
  filter(olink_suggested_signal == T) %>%
  mutate(untransformedNPX = (2^PCNormalizedNPX))

olink_5_w <- olink_5 %>%
  dplyr::select(SampleID, UniProt, PCNormalizedNPX) %>%
  pivot_wider(names_from = SampleID, values_from = PCNormalizedNPX, values_fill = NA) %>%
  dplyr::select(!contains("_r")) 
colnames(olink_5_w)[1] <- "PG.ProteinGroups"
colnames(olink_5_w) <- sub("O5-", "", colnames(olink_5_w))


## Make overlap list of protein groups for this experiment
f <- list(neat = neat_PG, 
          acid = acid_PG, 
          magnet = magnet_PG, 
          preomics = preomics_PG, 
          seer = seer_PG_NPs)

overlap <- list()

for (i in names(f)) {
  g <- f[[i]]
  h <- unique(g$PG.ProteinGroups)
  overlap[[i]] <- h
}
g <- olink_5_w 
h <- unique(g$PG.ProteinGroups)
overlap <- c(overlap, list(olink = h))

dd <- transpose(overlap, keep.names = "ID")

fwrite(dd, "data/processed/Exp5_UniProtIDs_ExpandedProteinGroups.csv", na = "", row.names = F)



#### Data Processing 

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


##* log2 transform ----

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


##* Impute (ImpSeq) ----

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


## t-testing
canc <- list(neat = neat_PG_imp, 
             acid = acid_PG_imp, 
             magnet = magnet_PG_imp, 
             preomics = preomics_PG_imp, 
             seer = seer_PG_imp,
             olink = olink_PG_imp)


volc_plot <- data.frame(PG.ProteinGroup = character(),
                        logfc = numeric(),
                        pvalue = numeric(),
                        qvalue = numeric(),
                        method = character())

for (i in names(canc)) {
  
  # choose a data frame
  all_df <- canc[[i]] %>%
    tibble::column_to_rownames(var = "PG.ProteinGroups") %>%
    dplyr::select(num_range("", 1:40))
  
  # make a grouping vector
  groups <- c(rep("Healthy", 20), rep("Cancer", 20))
  
  
  # run ROTS
  results <- ROTS(data = all_df, 
                  groups = groups, 
                  B = 5000, 
                  K = 2000, 
                  seed = 1234)
  
  summary(results, fdr = 0.05)
  
  # output df
  all_fc <- data.frame(results$logfc) %>%
    tibble::rownames_to_column(var = "PG.ProteinGroup") %>%
    mutate(logfc = results.logfc) %>%
    dplyr::select(-results.logfc) %>%
    mutate(pvalue = results$pvalue) %>%
    mutate(qvalue = results$FDR) %>%
    mutate(method = i)
  
  volc_plot <- rbind(volc_plot, all_fc)
  
}


volc_plot <- volc_plot %>%
  mutate(neglogpvalue = -log10(pvalue)) %>%
  mutate(diffexp = case_when(
    logfc >= 0.263 & qvalue <= 0.05 ~ "UP",
    logfc <= -0.263 & qvalue <= 0.05 ~ "DOWN",
    T ~ "NO"
  ))


#### Figure S16 ####
UP <- volc_plot %>%
  filter(diffexp == "UP")
DOWN <- volc_plot %>%
  filter(diffexp == "DOWN")

UPS <- UP %>%
  count(PG.ProteinGroup, name = "count") %>%
  arrange(desc(count))
DOWNS <- DOWN %>%
  count(PG.ProteinGroup, name = "count") %>%
  arrange(desc(count))
# proteins that are both up/down depending on method
BOTH <- UPS %>%
  inner_join(DOWNS, by = "PG.ProteinGroup")
BOTHE <- volc_plot %>%
  pivot_wider(
    id_cols = PG.ProteinGroup,
    names_from = method,
    values_from = diffexp
  ) %>%
  filter(PG.ProteinGroup %in% BOTH$PG.ProteinGroup)

UP_list <- UP %>%
  group_by(method) %>%
  summarise(ID_list = list(PG.ProteinGroup)) %>%
  deframe()
DOWN_list <- DOWN %>%
  group_by(method) %>%
  summarise(ID_list = list(PG.ProteinGroup)) %>%
  deframe()


pdf("reports/figures/Experiment5_UPDOWNoverlap_1,2FCCutoff_Alternative.pdf", width = 6, height = 3)
upset(fromList(UP_list), 
      sets = names(UP_list), 
      nintersects = NA,
      order.by = "freq",
      point.size = 1,
      line.size = 0.2,
      text.scale = 1,
      mainbar.y.label = 'Differential Proteins',
      sets.x.label = 'Differential Proteins per Method')

upset(fromList(DOWN_list), 
      sets = names(DOWN_list), 
      nintersects = NA,
      order.by = "freq",
      point.size = 1,
      line.size = 0.2,
      text.scale = 1,
      mainbar.y.label = 'Differential Proteins',
      sets.x.label = 'Differential Proteins per Method')
dev.off()