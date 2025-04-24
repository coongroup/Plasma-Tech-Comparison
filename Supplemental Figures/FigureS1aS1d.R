library(devtools)
library(tidyverse)
library(OlinkAnalyze)



pal <- c("#23B56D",
         "#F26439",
         "#5762AC",
         "#E38CBB",
         "#2BA7DF",
         "#99C13C")

#### Figure S1a ####
dire <- "data/raw/"
fi <- c("20240912_WFB_exp01_acid_3_0.csv",
        "20240912_WFB_exp01_magnet_3_0.csv",
        "20240912_WFB_exp01_neat_3_0.csv",
        "20240912_WFB_exp01_preomics_3_0.csv",
        "20240912_WFB_exp01_seer_3_NPB.csv",
        "20240912_WFB_exp01_seer_3_NPA.csv")

chromatograms <- NULL
for (i in fi) {
  
  #read csv
  chrom <- read.csv(paste0(dire, i))
  chrom <- tibble::rownames_to_column(chrom)
  chrom <- chrom[-c(1,2),]
  colnames(chrom) <- chrom[1,]
  chrom <- chrom[-1,]
  chrom$meth <- sapply(strsplit(i, "_"), function(x) paste(x[4], collapse = "_"))
  chrom$meth <- ifelse(grepl("NPA", i), "seerNPA",
                       ifelse(grepl("NPB", i), "seerNPB", chrom$meth))
  
  #make long data frame
  chromatograms <- rbind(chromatograms, chrom)
}

colnames(chromatograms) <- c("time", "abundance", "method")
chromatograms$time <- as.numeric(chromatograms$time)
chromatograms$abundance <- as.numeric(chromatograms$abundance)
chromatograms$method <- as.factor(chromatograms$method)
chromatograms$method <- factor(chromatograms$method, levels = c("neat", "acid", "preomics", "magnet", "seerNPA", "seerNPB"))

#### Plot Chromatorgrams 
labs <- chromatograms %>%
  dplyr::distinct(method) %>%
  dplyr::mutate(x = 0.5, y = 2.8)

#make plot
ggplot(chromatograms, aes(time, abundance/10000000000)) + 
  geom_line(size = 0.1) + 
  #scale_color_manual(values = col) +
  #ggtitle("Chromatograms") +
  xlab("Time (min)") +
  ylab("Abundance (x1e10)") +
  #ylim(-10, 0) +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 1, 2, 3)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 40, by = 8)) +
  guides(fill = guide_legend(title="Set")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank()) +
  facet_grid(method~ .) +
  geom_text(data = labs, aes(x = x, y = y, label = method), inherit.aes = FALSE, size = 2, color = "black", hjust = 0) 


#### Figure S1d ####
neat_PG <- read.csv("data/processed/PG_Matrix_Experiment1_neat.csv")
acid_PG <- read.csv("data/processed/PG_Matrix_Experiment1_acid.csv")
preomics_PG <- read.csv("data/processed/PG_Matrix_Experiment1_preomics.csv")
magnet_PG <- read.csv("data/processed/PG_Matrix_Experiment1_magnet.csv")
seer_PG <- read.csv("data/processed/PG_Matrix_Experiment1_seer.csv")

dop1 <- list(neat_PG = neat_PG, 
             acid_PG = acid_PG, 
             magnet_PG = magnet_PG, 
             preomics_PG = preomics_PG, 
             seer_PG = seer_PG)


for (i in names(dop1)) {
  
  dop1[[i]] <- dop1[[i]] %>% 
    select(-X, -R.Condition) %>%
    pivot_longer(cols = -PG.ProteinGroups, names_to = "group", values_to = "value") %>%
    filter(grepl("_1|_2|_3|_4|_5", group)) %>%
    filter(!grepl("_0_2|_1_2", group)) %>%
    filter(!grepl("_2024", group)) %>%
    pivot_wider(names_from = group, values_from = value)
  
}

neat_PG <- dop1[["neat_PG"]]
acid_PG <- dop1[["acid_PG"]]
magnet_PG <- dop1[["magnet_PG"]]
preomics_PG <- dop1[["preomics_PG"]]
seer_PG <- dop1[["seer_PG"]]



#split NPA and NPB
NPA_PG <- seer_PG[, c("PG.ProteinGroups", grep("NPA", names(seer_PG), value = TRUE))]
NPB_PG <- seer_PG[, c("PG.ProteinGroups", grep("NPB", names(seer_PG), value = TRUE))]

#combine seer NPA and NPB
#use this to take either NPA or NPB for all samples for each pg, depending on which set is more complete
#make a filter to choose whichever row has fewer missing values for NPA or NPB. If same, choose NPA.
NPA <- apply(seer_PG[,grep("NPA", names(seer_PG))], 1, function(x) table(unlist(unname(x)) > 0)[1])
NPA[is.na(NPA)] <-0

NPB <- apply(seer_PG[,grep("NPB", names(seer_PG))], 1, function(x) table(unlist(unname(x)) > 0)[1])
NPB[is.na(NPB)] <-0

filter_NP <- NPA>=NPB
table(is.na(filter_NP))

#separate NPA and NPB into different dataframes
NPA_df <- seer_PG[,c(1,grep("NPA", names(seer_PG)))]
NPB_df <- seer_PG[,c(1,grep("NPB", names(seer_PG)))]

#apply filter to choose correct rows for either dataframe
NPA_df_filter <- NPA_df[filter_NP,]
NPB_df_filter <- NPB_df[!filter_NP,]
names(NPA_df_filter) <- sub(pattern = "_NPA", replacement = "", x = names(NPA_df_filter))
names(NPB_df_filter) <- sub(pattern = "_NPB", replacement = "", x = names(NPB_df_filter))
names(NPA_df_filter) <- sub(pattern = "_20240915111008", replacement = "_2", x = names(NPA_df_filter))
names(NPB_df_filter) <- sub(pattern = "_20240915181532", replacement = "_2", x = names(NPB_df_filter))
NPA_df_filter_1 <- NPA_df_filter[,c(1, order(names(NPA_df_filter)[-1])+1)]
NPB_df_filter_1 <- NPB_df_filter[,c(1, order(names(NPB_df_filter)[-1])+1)]

seer_PG_NPs <- rbind(NPA_df_filter, NPB_df_filter)
#save combined quant matrix
write.csv(seer_PG_NPs, file = "data/processed/PG_Matrix_Experiment1_NPs_separate.csv")
#need to add this one manually


dop <- list(neat_PG = neat_PG, 
            acid_PG = acid_PG, 
            magnet_PG = magnet_PG, 
            preomics_PG = preomics_PG, 
            NPA_PG = NPA_PG,
            NPB_PG = NPB_PG)

metrics <- data.frame(
  Method = character(),
  Total = numeric(),
  Complete = numeric(),
  Average = numeric(),
  StDev = numeric()
)

for (i in names(dop)) {
  
  x <- dop[[i]]
  
  #remove the first run and R.Condition from each of the data frames
  x <- x[,-c(1, 3:4)]
  
  #make sure there is at least one value in the measurement columns
  x <- x[rowSums(is.na(x[, 2:6])) != 5, ]
  
  #all PGs
  n <- nrow(x) 
  
  #complete cases
  cc <- nrow(na.omit(x))
  
  #average + sd
  mean <- mean(colSums(!is.na(x[, 2:6]))) 
  sd <- sd(colSums(!is.na(x[, 2:6]))) 
  
  temp <- data.frame(Method = i,
                     Total = n,
                     Complete = cc,
                     Average = mean,
                     StDev = sd)
  metrics <- rbind(metrics, temp)
}

#add the olink metrics by rbind
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
ol <- data.frame(Method = "Olink_PG",
                 Total = nrow(olink_w),
                 Complete = nrow(na.omit(olink_w)),
                 Average = mean(colSums(!is.na(olink_w[, 2:6]))),
                 StDev = sd(colSums(!is.na(olink_w[, 2:6]))))
metrics <- rbind(metrics, ol)

#add the seer NP combined by rbind
seer_PG_NPs <- seer_PG_NPs[,-3]
seer_PG_NPs <- seer_PG_NPs[rowSums(is.na(seer_PG_NPs[, 2:6])) != 5, ]
n <- nrow(seer_PG_NPs) 
cc <- nrow(na.omit(seer_PG_NPs))
mean <- mean(colSums(!is.na(seer_PG_NPs[, 2:6]))) 
sd <- sd(colSums(!is.na(seer_PG_NPs[, 2:6]))) 
temp <- data.frame(Method = "seer_PG_NPs",
                   Total = n,
                   Complete = cc,
                   Average = mean,
                   StDev = sd)
metrics <- rbind(metrics, temp)


metrics_s <- metrics %>%
  filter(Method %in% c("NPA_PG", "NPB_PG", "seer_PG_NPs"))

metrics_s_long <- metrics_s %>%
  pivot_longer(
    cols = c(Average, Total, Complete),     
    names_to = "Metric",          
    values_to = "Value" 
  ) %>%
  mutate(StDev = ifelse(Metric == "Average", Value, NA))

metrics_s_long$Metric <- factor(metrics_s_long$Metric, levels = c("Complete", "Average", "Total"))

#plot bars 
ggplot(metrics_s_long, aes(Method, Value, alpha = Metric)) + 
  geom_bar(position = position_dodge(width = 0.8), 
           stat = "identity",
           width = 0.7,
           fill = pal[5]) + 
  geom_text(label = round(metrics_s_long$Value, 0), 
            angle = 0, 
            color = "black",
            vjust = -0.5, 
            hjust = 0.5,
            size = 2, 
            position = position_dodge(width = 0.8)) +
  labs(x = NULL,
       y = "Protein Groups") +
  scale_alpha_manual(values = c(1, 0.75, 0.5)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = "bottom", 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm")
  ) 

