
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)

#Colors#
#Make a classic palette
col <- brewer.pal(8, "Set2") 



spec_PG <- read.csv("data/processed/PG_Matrix_Experiment1_separate_combined_noNPANPB.csv")

spec_PG_PCA <- spec_PG %>% 
  pivot_longer(cols = -PG.ProteinGroups, names_to = "group", values_to = "value") %>%
  filter(grepl("_1|_2|_3|_4|_5", group)) %>%
  filter(!grepl("_0_2|_1_2", group)) %>%
  pivot_wider(names_from = group, values_from = value)

spec_PG_PCA <- merge(spec_PG_PCA, olink_w, by.x = "PG.ProteinGroups", by.y = "UniProt", all = T)

selec <- c("neat", "acid", "preomics", "magnet", "seer", "O1")

t_variance <- data.frame(prop = numeric(),
                         method = character())


t_scores <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_scores) <- paste0("PC", 1:5)
t_scores[] <- lapply(t_scores, as.numeric)
t_scores["method"] <- lapply("method", function(x) character(0))
t_scores["Samples"] <- lapply("Samples", function(x) character(0))


for (i in selec) {
  
  # choose a data frame
  pca_set <- spec_PG_PCA %>%
    select(PG.ProteinGroups, contains(i)) %>%
    filter(complete.cases(.))
  
  # transpose
  t_pca_set <- as.data.frame(t(pca_set))
  
  # get them into shape
  colnames(t_pca_set) <- as.character(t_pca_set[1, ])
  t_pca_set <- t_pca_set[-1, ]
  t_pca_set_n <- rownames(t_pca_set)
  t_pca_set_1 <- as.data.frame(lapply(t_pca_set, function(x) as.numeric(as.character(x))))
  t_pca_set_1 <- log2(t_pca_set_1)
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
    mutate(Samples = t_pca_set$Samples) %>%
    mutate(method = i)
  
  t_scores <- rbind(t_scores, scores)
  
}

t_scores <- t_scores %>%
  mutate(order = sapply(strsplit(Samples, "_"), function(x) paste(x[5], collapse = "_"))) %>%
  mutate(order = if_else(order == "NA", substr(Samples, 4, 4), order))


#### Figure S6a ####
ggplot(t_scores, aes(PC1, PC2)) + 
  geom_point(shape = 16,
             size = 1,
             aes(color = order)) +
  scale_color_viridis(option = "plasma", discrete = T) +
  xlab("PC1") +
  ylab("PC2") +
  xlim(c(-80, 80)) +
  ylim(c(-80, 80)) +
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
ggsave("reports/figures/Experiment1_All_CompleteCases_order_PCA_small.pdf", 
       width = 16, height = 10, units = "cm")


#### Figure S6b ####
t_variance <- t_variance %>%
  mutate(pro = prop*100) %>%
  group_by(method) %>%
  mutate(index = row_number()) %>%
  mutate(method = factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer", "O1")))


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
ggsave("reports/figures/Experiment1_All_CompleteCases_ScreePlot_small.pdf", 
       width = 8.8, height = 4, units = "cm")



#### Figure S6c ####
# using dataframes from Exp5_HealthyvsCancer script
heats <- spec_PG_PCA %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "colnames", values_to = "values") %>%
  mutate(sample = word(colnames, 4, 5, sep = "_")) %>%
  mutate(sample = case_when(
    colnames == "O1-1" ~ "olink_1",
    colnames == "O1-2" ~ "olink_2",
    colnames == "O1-3" ~ "olink_3",
    colnames == "O1-4" ~ "olink_4",
    colnames == "O1-5" ~ "olink_5",
    T ~ sample
  )) %>%
  select(-colnames) %>%
  pivot_wider(names_from = "sample", values_from = "values")




selec <- c("neat", "acid", "preomics", "magnet", "seer", "olink")

pheatmaps <- list()
for (i in selec) {
  
  # choose a data frame
  ne <- heats %>%
    select(PG.ProteinGroups, contains(i)) %>%
    mutate(across(2:6, log2)) %>%
    filter(rowSums(!is.na(.[,2:6])) > 0)
  
  nea <- ne %>%
    mutate(across(2:6, ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) # column mean
  #    mutate(across(2:6, ~ ifelse(is.na(.), rowMeans(ne[, 2:6], na.rm = TRUE), .))) #row mean
  
  pheats <- 
    pheatmap(ne[,2:6],
             clustering_distance_rows = dist(nea[,2:6]), 
             clustering_distance_cols = dist(t(nea[,2:6])),
             #color = viridis(24, direction = 1, option = "plasma"),
             color = rev(colorRampPalette(brewer.pal(24, "RdYlBu"))(24)),
             na_col = "white",
             cluster_rows = T,
             cluster_cols = T,
             treeheight_row = 0,
             treeheight_col = 10,
             show_rownames = F,
             show_colnames = T,
             fontsize = 16,
             border_color = NA,
             main = paste(i, "unimputed"),
             width = 4,
             height = 3)
  
  pheatmaps[[paste0("heatmap_unimputed_", i)]] <- pheats$gtable
  
}

pdf("reports/figures/Experiment1_All_Unimputed_Heatmaps_columnmean.pdf", width = 24, height = 8)
grid.arrange(
  grobs = pheatmaps,
  nrow = 1)
dev.off()

png("reports/figures/Experiment1_All_Unimputed_Heatmaps_columnmean.png", width = 2400, height = 800)
grid.arrange(
  grobs = pheatmaps,
  nrow = 1)
dev.off()
