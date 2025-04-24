#### libraries/colors/files ####
library(devtools)
library(tidyverse)
library(rrcovNA)
library(ROTS)



# files #

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



#### Figure 4a ####

canc <- list(neat = neat_PG, 
             acid = acid_PG, 
             magnet = magnet_PG, 
             preomics = preomics_PG, 
             seer = seer_PG_NPs,
             olink = olink_5_w)

miss <- data.frame(
  Method = character(),
  Order = integer(),
  Value = numeric()
)

lines <- data.frame(
  Method = character(),
  All = integer(),
  Half = integer(),
  Min = integer()
)


for (i in names(canc)) {
  
  z <- canc[[i]]
  
  #Calc percent of samples each protein was detected in
  percent_non_na <- rowMeans(!is.na(z)) * 100
  non_na_counts_ordered <- percent_non_na[order(-percent_non_na)]
  
  all <- sum(percent_non_na == 100) 
  half <- sum(percent_non_na >= 50)
  min <- sum(percent_non_na >= 0)
  
  #make missing plotting dataframe
  missing <- data.frame(Method = i,
                        Order = seq_along(non_na_counts_ordered), 
                        Value = non_na_counts_ordered)
  
  #make vertical line dataframe
  ls <- data.frame(Method = i,
                   All = all,
                   Half = half,
                   Min = min)
  
  
  miss <- rbind(miss, missing)
  lines <- rbind(lines, ls)
  
}

labse <- miss %>%
  dplyr::distinct(Method) %>%
  dplyr::mutate(x = 6000, y = 90)

miss$Method <- factor(miss$Method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink"))

ggplot(miss, aes(Order, Value)) + 
  geom_line(size = 0.2,
            color = "black",
            fill = col[8]) + 
  #geom_point(size = 1,
  #           color = "black",
  #           fill = NA,
  #          shape = 21) +
  geom_segment(data = lines, 
               aes(x = All,
                   xend = All,
                   y = 0,
                   yend = 100), 
               linetype = "solid", 
               color = "gray", 
               linewidth = 0.2) +
  geom_segment(data = lines, 
               aes(x = Half,
                   xend = Half,
                   y = 0,
                   yend = 50), 
               linetype = "solid", 
               color = "gray", 
               linewidth = 0.2) +
  #scale_fill_manual(values = col[5]) +
  xlab("Protein Groups") +
  ylab("% of Samples Quantified") +
  #ylim(0, 20) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,1)) +
  #coord_cartesian(xlim = c(0, 10000), ylim = c(0, 101)) +
  #guides(fill = guide_legend(title="Impairment\nStatus")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(1, "lines"), 
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 12),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank()) +
  facet_grid(factor(Method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink")) ~ ., scales = "free") +
  geom_text(data = lines, 
            aes(x = Min, y = 10, label = Min), 
            inherit.aes = FALSE, 
            size = 3, 
            color = "black", 
            hjust = 3) +
  geom_text(data = lines, 
            aes(x = Half, y = 60, label = Half), 
            inherit.aes = FALSE, 
            size = 3, 
            hjust = 0) +
  geom_text(data = lines, 
            aes(x = All, y = 90, label = All), 
            inherit.aes = FALSE, 
            size = 3, 
            hjust = 0) +
  geom_text(data = labse, 
            aes(x = x, y = y, label = Method),
            size = 3, 
            inherit.aes = FALSE)








#### Figure 4b ####

##* adjust column names ----
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


##* Olink Data Median Normalization ----
olink_PG_imp <- olink_PG_imp %>%
  mutate(across(where(is.numeric), ~ . - median(., na.rm = TRUE)))


##* ROTS DEA ----

##* ttesting and plotting ----

#switch with imputed and non imputed values here
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

counts <- volc_plot %>%
  group_by(method) %>%
  filter(diffexp %in% c("DOWN", "UP", "NO")) %>%
  count(diffexp)

volc_plot$method <- factor(volc_plot$method,
                           levels = c("neat", "acid", "preomics", "magnet", "seer", "olink"))


ggplot(volc_plot, aes(logfc, neglogpvalue, color = diffexp, size = diffexp)) + 
  geom_point() +
  geom_vline(xintercept=c(-0.263, 0.263), 
             col="black",
             size = 0.2) +
  #geom_hline(yintercept=-log10(0.05), 
  #           col="black",
  #           size = 0.2) +
  scale_color_manual(values = c("#9e1b45", col[8], "#9e1b45")) +
  scale_size_manual(values = c(0.5,0.1,0.5)) +
  scale_x_continuous(limits = c(-max(abs(volc_plot$logfc)), max(abs(volc_plot$logfc))), breaks = seq(-5, 5, by = 2)) +
  #geom_text_repel(data = subset(volc_plot, diffexp != "NO"), aes(label = gene), size = 2) +
  xlab("Log2 Fold Change (Cancer/Healthy)") +
  ylab("-Log10 Adjusted P-Value") +
  #xlim(-2, 2) +
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
        legend.position = "bottom") +
  facet_wrap(. ~ factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink")),
             scales = "free",
             ncol = 2) +
  geom_text(data = volc_plot %>% distinct(method),
            aes(x = -3.5, y = 8, label = method),
            size = 3, 
            color = "black", 
            hjust = 0) +
  geom_text(data = counts[counts$diffexp == "DOWN",], 
            aes(x = -3, y = Inf, 
                label = paste(n)),
            hjust = 1.1, vjust = 1.5, size = 3, show.legend = FALSE) +
  geom_text(data = counts[counts$diffexp == "UP",], 
            aes(x = 3, y = Inf, 
                label = paste(n)),
            hjust = 1.1, vjust = 1.5, size = 3, show.legend = FALSE) +
  geom_text(data = counts[counts$diffexp == "NO",], 
            aes(x = 0, y = Inf, 
                label = paste(n)),
            hjust = 0.5, vjust = 1.5, size = 3, show.legend = FALSE) 

