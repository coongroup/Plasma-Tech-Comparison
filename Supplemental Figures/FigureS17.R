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


#### Figure S17 ####
library(enrichR)
setEnrichrServer("https://maayanlab.cloud/Speedrichr/")
library(UniProt.ws)
up <- UniProt.ws(taxId = 9606)

dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
         "GO_Biological_Process_2023")

#subset data by fold change and pvalue, and add gsea to a long data frame with every method

methodss <- c("neat", "acid", "preomics", "magnet", "seer", "olink")

enrichr <- data.frame(Term = character(),
                      Rank = character(),
                      P.Value = numeric(),
                      Adjusted.P.Value = numeric(),
                      Old.P.Value = integer(),
                      Old.Adjusted.P.Value = integer(),
                      Odds.Ratio = numeric(),
                      Combined.Score = numeric(),
                      Genes = character(),
                      method = character())

#split protein groups and olink groups into individual rows
volc_plot1 <- volc_plot %>%
  separate_rows(PG.ProteinGroup, sep = ";") %>%
  separate_rows(PG.ProteinGroup, sep = "_")

for (i in methodss) {
  
  gene_list <- volc_plot1 %>%
    filter(method == i) 
  mapping <- UniProt.ws::select(up, keys = gene_list$PG.ProteinGroup, keytype = "UniProtKB", columns = "gene_primary")
  gene_list <- gene_list %>%
    left_join(mapping, by = c("PG.ProteinGroup" = "From")) 
  
  background <- gene_list$Gene.Names..primary.
  background <- unlist(strsplit(background, ";"))
  write.csv(background, "data/processed/Exp5_EnrichR_neat_background.csv")
  
  gene_list <- gene_list %>%
    filter(diffexp != "NO")
  
  input <- gene_list$Gene.Names..primary.
  input <- unlist(strsplit(input, ";"))
  write.csv(input, "data/processed/Exp5_EnrichR_neat_input.csv")
  
  enriched1 <- enrichr(input, dbs, background = background)
  
  asd <- enriched1$GO_Biological_Process_2023 %>%
    mutate(method = i)
  
  enrichr <- rbind(enrichr, asd)
  
}

enrichr <- enrichr %>%
  mutate(`-log10qvalue` = -log10(Adjusted.P.value)) %>%
  mutate(Rank = as.numeric(enrichr$Rank))

ggplot(enrichr, aes(Odds.Ratio, `-log10qvalue`)) +
  geom_point(shape = 16,
             size = 2,
             aes(color = Rank)) +
  scale_color_viridis(option = "plasma") +
  xlab("Odds Ratio") +
  ylab("-log10(q-value)") +
  #xlim(c(-150, 150)) +
  ylim(c(0, 30)) +
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