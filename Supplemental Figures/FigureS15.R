
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



#### Figure S15 ####


canc <- list(neat = neat_PG_imp, 
             acid = acid_PG_imp, 
             magnet = magnet_PG_imp, 
             preomics = preomics_PG_imp, 
             seer = seer_PG_imp,
             olink = olink_PG_imp)

t_variance <- data.frame(prop = numeric(),
                         method = character())


t_scores <- data.frame(matrix(ncol = 40, nrow = 0))
colnames(t_scores) <- paste0("PC", 1:40)
t_scores[] <- lapply(t_scores, as.numeric)
t_scores["Samples"] <- lapply("Samples", function(x) numeric(0))
t_scores[c("group", "method")] <- lapply(c("group", "method"), function(x) character(0))

for (i in names(canc)) {
  
  # choose a data frame
  pca_set <- canc[[i]]
  
  # transpose
  t_pca_set <- as.data.frame(t(pca_set))
  
  # get them into shape
  colnames(t_pca_set) <- as.character(t_pca_set[1, ])
  t_pca_set <- t_pca_set[-1, ]
  t_pca_set_n <- rownames(t_pca_set)
  t_pca_set_1 <- as.data.frame(lapply(t_pca_set, function(x) as.numeric(as.character(x))))
  rownames(t_pca_set_1) <- t_pca_set_n
  t_pca_set <- t_pca_set_1
  t_pca_set <- tibble::rownames_to_column(t_pca_set, "Samples")
  
  # PCA
  pca_score <- prcomp(t_pca_set[,c(2:ncol(t_pca_set))],
                      scale. = T)
  
  # make variables to plot
  #scree
  explained_variance <- pca_score$sdev^2 / sum(pca_score$sdev^2)
  
  variance <- data.frame(prop = explained_variance,
                         method = i)
  
  t_variance <- rbind(t_variance, variance)
  
  #pca scores
  scores <- as.data.frame(pca_score$x)
  scores <- scores %>%
    mutate(Samples = as.numeric(t_pca_set$Samples)) %>%
    mutate(group = case_when(
      Samples >= 1 & Samples <= 20 ~ "healthy",   
      Samples >= 21 & Samples <= 40 ~ "cancer"),
      method = i)
  
  t_scores <- rbind(t_scores, scores)
  
}


# PCA Plots
ggplot(t_scores, aes(PC1, PC2)) + 
  geom_point(shape = 16,
             size = 1,
             aes(color = group)) +
  scale_color_manual(values = c("#9f1c46", col[7])) +
  xlab("PC1") +
  ylab("PC2") +
  xlim(c(-150, 150)) +
  ylim(c(-100, 100)) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0),
        axis.ticks = element_line(size = 0.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        legend.position = "bottom",
        legend.justification = "left",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(color = "gray", fill = NA, size = 0.2)) +
  facet_wrap(factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink")) ~ .,
             scales = "free") 
ggsave("reports/figures/Experiment5_All_50pMissingImputed_olinknormalized_PCA_small.pdf", 
       width = 18, height = 12, units = "cm")

# scree plot
t_variance <- t_variance %>%
  mutate(pro = prop*100) %>%
  group_by(method) %>%
  mutate(index = row_number()) %>%
  mutate(method = factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink")))


ggplot(t_variance, aes(index, pro), fill = method, color = method) +
  geom_point(shape = 16,
             size = 1,
             aes(color = method)) +
  geom_line(aes(color = method),
            alpha = 0.5) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlab("Component") +
  ylab("Proportion of Variance") +
  xlim(1, 5) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(1, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")) 