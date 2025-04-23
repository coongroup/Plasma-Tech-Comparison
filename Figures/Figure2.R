
library(devtools)
library(tidyverse)
library(OlinkAnalyze)



pal <- c("#23B56D",
         "#F26439",
         "#5762AC",
         "#E38CBB",
         "#2BA7DF",
         "#99C13C")



#### files ####
spec_PG <- read.csv("data/processed/PG_Matrix_Experiment1_separate_combined_noNPANPB.csv")
file_info <- read.csv("data/metadata/Experiment1_file_info.csv")
file_info <- file_info[-1]
file_info$run <- sub(pattern = ".raw", replacement = "", x = file_info$Runs)

#make function to calculate rsd
rsd <- function(featureRow, index = 1:length(featureRow)){sd(featureRow[index], na.rm = T)/mean(featureRow[index], na.rm = T) * 100}


## REMOVE Replicate run FOR EACH METHOD
spec_PG <- spec_PG %>% 
  pivot_longer(cols = -PG.ProteinGroups, names_to = "group", values_to = "value") %>%
  filter(grepl("_1|_2|_3|_4|_5", group)) %>%
  filter(!grepl("_0_2|_1_2", group)) %>%
  pivot_wider(names_from = group, values_from = value)


## Open Olink File 
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

#add column of RSD across all QCs
olink_w_rsd <- cbind(olink_w, rsd = apply(olink_w[,grepl("O1",colnames(olink_w))], 1, rsd))
olink_w_rsd$Method <- "Olink"



## CVs Analysis 
#add column of RSD across all QCs
spec_PG_rsd <- cbind(spec_PG, rsd_all = apply(spec_PG[,grepl("WFB",colnames(spec_PG))], 1, rsd))

#add columns of RSDs within plates
spec_PG_rsd$rsd_neat <- apply(spec_PG_rsd[, grepl("neat",colnames(spec_PG_rsd))], 1, rsd)
spec_PG_rsd$rsd_acid <- apply(spec_PG_rsd[, grepl("acid",colnames(spec_PG_rsd))], 1, rsd)
spec_PG_rsd$rsd_magnet <- apply(spec_PG_rsd[, grepl("magnet",colnames(spec_PG_rsd))], 1, rsd)
spec_PG_rsd$rsd_preomics <- apply(spec_PG_rsd[, grepl("preomics",colnames(spec_PG_rsd))], 1, rsd)
spec_PG_rsd$rsd_seer <- apply(spec_PG_rsd[, grepl("seer",colnames(spec_PG_rsd))], 1, rsd)

spec_PG_rsd <- merge(spec_PG_rsd, olink_w_rsd, by.x = "PG.ProteinGroups", by.y = "UniProt", all = T)



#### Figure 2a-b ####

#limit to complete cases for plotting
#spec_PG_rsd <- spec_PG_rsd[complete.cases(spec_PG_rsd),]

#Make a long dataframe for plotting # use this for the next 2 plots
spec_PG_rsd_long <- spec_PG_rsd[,-39] %>%
  pivot_longer(cols = -PG.ProteinGroups, names_to = "Group", values_to = "Value")
spec_PG_rsd_long$Group <- sub(pattern = "X", replacement = "", x = spec_PG_rsd_long$Group)


#boxplot of RSD for each run in run order
spec_PG_rsd_long <- spec_PG_rsd_long[grepl("rsd", spec_PG_rsd_long$Group),]
spec_PG_rsd_long <- spec_PG_rsd_long[!grepl("all", spec_PG_rsd_long$Group),]
spec_PG_rsd_long <- spec_PG_rsd_long %>%
  filter(!is.na(Value))

spec_PG_rsd_long$Group <- sub("rsd_", "", spec_PG_rsd_long$Group)
spec_PG_rsd_long$Group <- sub("rsd", "Olink", spec_PG_rsd_long$Group)


#boxplot of RSD for each run in run order
spec_PG_rsd_long <- spec_PG_rsd_long[grepl("rsd", spec_PG_rsd_long$Group),]
spec_PG_rsd_long <- spec_PG_rsd_long[!grepl("all", spec_PG_rsd_long$Group),]
spec_PG_rsd_long <- spec_PG_rsd_long %>%
  filter(!is.na(Value))


medians <- spec_PG_rsd_long %>%
  group_by(Group) %>%
  summarise(median_value = median(Value, na.rm = T))

spec_PG_rsd_long$Group <- factor(spec_PG_rsd_long$Group, 
                                 levels = c("rsd_neat", "rsd_acid", "rsd_preomics", "rsd_magnet", "rsd_seer", "rsd"))

ggplot(spec_PG_rsd_long, aes(Group, Value, fill = Group)) + 
  geom_violin(alpha = 0.5, width = 1) +
  geom_boxplot(width = 0.08, alpha = 0.5) +
  #scale_fill_manual(values = mycolors) +
  scale_fill_manual(values = col) +
  geom_text(data = medians, aes(x = Group, y = median_value, label = round(median_value, 2)),
            vjust = -40, color = "black", size = 4) +
  ggtitle("%CV PG") +
  #xlab("Run Order") +
  ylab("CV") +
  #ylim(4000, 7000) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,50,100,150,200)) +
  #scale_x_continuous(expand = c(0,1)) +
  guides(fill = guide_legend(title="Plate")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = 'bold', hjust = 0.5), 
        axis.title = element_text(size = 24, face = 'bold'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12)) 
ggsave("reports/figures/Experiment1_PG_CVBoxplots_separate_combined_olink_CompleteCases1_noNPANPB.pdf", width = 32, height = 16, units = "cm")



#### Figure 2c ####
share <- spec_PG_rsd[complete.cases(spec_PG_rsd),][1]

meth <- c("neat", "acid", "preomics", "magnet", "seer")
quants <- data.frame(method = character(),
                     ProteinGroup = character(),
                     Median = numeric(),
                     rsd = numeric())
for (i in meth) {
  #make a name
  n <- paste0(i, "_median")
  
  medians <- apply(spec_PG[, grepl(i,colnames(spec_PG))], 1, median, na.rm = TRUE)
  
  #make dataframe
  quant <- cbind(i, spec_PG$PG.ProteinGroups, medians, spec_PG_rsd[[paste0("rsd_",i)]])
  quant <- as.data.frame(quant)
  quant$medians <- as.numeric(quant$medians)
  quant$V4 <- as.numeric(quant$V4)
  colnames(quant) <- c("method", "ProteinGroup", "Median", "rsd")
  
  quants <- rbind(quants, quant)
  
}
quants$shared <- ifelse(quants$ProteinGroup %in% share$PG.ProteinGroups, "yes", "no")

lab <- quants %>%
  dplyr::distinct(method) %>%
  dplyr::mutate(x = 0, y = 230)
poin <- quants %>% 
  group_by(method) %>%
  summarize(non_na_count = sum(!is.na(Median)))

#make plot
ggplot(quants, aes(log2(Median), rsd, fill = shared, alpha = shared)) + 
  geom_point(shape = 21, 
             size = 2, 
             stroke = 0) + 
  #geom_smooth( se = FALSE, color = "red") +
  #scale_fill_manual(values = mycolors) +
  scale_alpha_manual(values = c(0.1, 0.2)) +
  scale_fill_manual(values = c("black", col[4])) +
  xlab("log2(Median PG Quantity)") +
  ylab("CV") +
  #ylim(4000, 7000) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 240)) +
  #scale_x_continuous(expand = c(0,1)) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(1, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank()) +
  facet_wrap( . ~ factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer")),
              scales = "free") +
  geom_text(data = lab, 
            aes(x = x, y = y, label = method), 
            inherit.aes = FALSE, 
            size = 2, 
            color = "black", 
            hjust = 0) +
  geom_text(data = poin, 
            aes(x = 21, y = 230, label = paste("n =", non_na_count)), 
            inherit.aes = FALSE, 
            size = 2, 
            hjust = 0) 
ggsave(paste0("reports/figures/Experiment1_PGQuantvsRSD_noNPANPB_small.pdf"), width = 7.5, height = 4, units = "in")


#separate graph for Olink
#make wider to calculate RSD
olink_w <- olink_1 %>%
  select(SampleID, UniProt, untransformedNPX) %>%
  pivot_wider(names_from = SampleID, values_from = untransformedNPX, values_fill = NA) %>%
  select(UniProt, `O1-1`, `O1-2`, `O1-3`, `O1-4`, `O1-5`)

#add column of RSD across all QCs
olink_w_rsd <- cbind(olink_w, rsd = apply(olink_w[,grepl("O1",colnames(olink_w))], 1, rsd))
olink_w_rsd$Method <- "Olink"

medi <- apply(olink_w[, grepl("O1",colnames(olink_w))], 1, median, na.rm = TRUE)

#make dataframe
quan <- cbind(olink_w_rsd$Method, 
              olink_w$UniProt, 
              medi, 
              olink_w_rsd$rsd)
quan <- as.data.frame(quan)
quan$medi <- as.numeric(quan$medi)
quan$V4 <- as.numeric(quan$V4)
colnames(quan) <- c("method", "ProteinGroup", "Median", "rsd")

quan$shared <- ifelse(quan$ProteinGroup %in% share$PG.ProteinGroups, "yes", "no")

#add HPPP abundance for the different Olink proteins for plotting
HPPP <- read.csv("data/metadata/2021_HPPP_Protein_List.csv") %>% mutate(ProteinGroup = protein_identifier) %>% select(ProteinGroup, plasma_abundance)

quan <- quan %>%
  left_join(HPPP, by = "ProteinGroup")

labe <- quan %>%
  dplyr::distinct(method) %>%
  dplyr::mutate(x = 0, y = 220)
poine <- quan %>% 
  group_by(method) %>%
  summarize(non_na_count = sum(!is.na(plasma_abundance)))


#make plot
ggplot(quan, aes(plasma_abundance, rsd, fill = shared, alpha = shared)) + 
  geom_point(shape = 21, 
             size = 2, 
             stroke = 0) + 
  #geom_smooth( se = FALSE, color = "red") +
  #scale_fill_manual(values = mycolors) +
  scale_alpha_manual(values = c(0.1, 0.2)) +
  scale_fill_manual(values = c("black", col[4])) +
  xlab("HPPP abundance") +
  ylab("CV") +
  #ylim(4000, 7000) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 240)) +
  #scale_x_continuous(expand = c(0,1)) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.position = "none",
        axis.text.x = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks.x = element_line(size = 0.2),
        axis.ticks.y = element_blank(),
        strip.text = element_blank()) +
  geom_text(data = labe, 
            aes(x = x, y = y, label = method), 
            inherit.aes = FALSE, 
            size = 2, 
            color = "black", 
            hjust = 0) +
  geom_text(data = poine, 
            aes(x = 5, y = 220, label = paste("n =", non_na_count)), 
            inherit.aes = FALSE, 
            size = 2, 
            hjust = 0) 
ggsave(paste0("reports/figures/Experiment1_PGQuantvsRSD_noNPANPB_olink_small.pdf"), width = 2.1, height = 2.1, units = "in")




