library(devtools)
library(tidyverse)
library(OlinkAnalyze)



pal <- c("#23B56D",
         "#F26439",
         "#5762AC",
         "#E38CBB",
         "#2BA7DF",
         "#99C13C")



#### Figure 3a ####
# read in default report
exp2_PG <- read.csv("data/processed/PG_Matrix_Experiment2_separate_combined_noNPANPB.csv")
exp2_pep <- read.csv("data/processed/pep_Matrix_Experiment2_separate_combined.csv")
all <- read.delim("data/spectronaut_output/20240928_223855_WFB_Report (Normal)_Experiment2_CombinedSNE.tsv", sep = "\t")

# CRP
all <- all %>%
  filter(PG.ProteinGroups == "P02741")
crp <- unique(all$EG.ModifiedSequence)

#limit tables to CRP PGs or Peps
PG_df_CRP <- exp2_PG %>%
  filter(PG.ProteinGroups == "P02741")

#make longer
PG_df_CRP_all_long <- PG_df_CRP[, -c(1)] %>%
  pivot_longer(cols = everything(), names_to = "Run", values_to = "Value")

#add columns for grouping to do average and sd
PG_df_CRP_all_long$rep <- sapply(strsplit(PG_df_CRP_all_long$Run, "_"), function(x) paste(x[4:6], collapse = "_"))
PG_df_CRP_all_long$com <- gsub("[0-9]+", "", PG_df_CRP_all_long$rep)
#log2 transform
PG_df_CRP_all_long$log <- log2(PG_df_CRP_all_long$Value)

#take mean and sd of the groups of 3
summary_PG_CRP <- PG_df_CRP_all_long %>%
  group_by(com) %>%
  summarise(
    mean_value = mean(log),
    sd_value = sd(log),
    se_value = sd_value / sqrt(n())  # Standard Error
  )

#add a column with the CRP dilutions depending on letter
summary_PG_CRP$dil <- ifelse(grepl("a_", summary_PG_CRP$com), 1,
                             ifelse(grepl("b", summary_PG_CRP$com), 2,
                                    ifelse(grepl("_c", summary_PG_CRP$com), 5,
                                           ifelse(grepl("_d", summary_PG_CRP$com), 10,
                                                  ifelse(grepl("_e", summary_PG_CRP$com), 100, NA)))))

summary_PG_CRP$method <- sapply(strsplit(summary_PG_CRP$com, "_"), function(x) paste(x[1], collapse = "_"))
summary_PG_CRP$method <- factor(summary_PG_CRP$method, levels = c("neat", "acid", "preomics", "magnet", "seer"))

#make fit data
fit <- lm(mean_value ~ log2(dil), data = summary_PG_CRP)


ggplot(summary_PG_CRP, aes(log2(dil), mean_value)) + 
  geom_point(shape = 21, 
             size = 1, 
             stroke = 0,
             fill = col[3]) + 
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                width = 0.2,
                size = 0.2) +
  geom_smooth(method = "lm", se = F, color = "gray",
              size = 0.2,
              linetype = "dashed") +
  stat_poly_eq(
    aes(label = ..rr.label..), 
    size = 2,
    formula = y ~ x,  
    label.x = "right",  
    label.y = "top",
    rr.digits = 4
  ) +
  xlab("log2(Dilution) 1x, 2x, 5x, 10x, 100x") +
  ylab("log2(CRP Quant)") +
  ylim(5, 24) +
  #scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(expand = c(0,1)) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        strip.text = element_blank()) +
  facet_grid(. ~ factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer")),
             scales = "free") +
  geom_text(data = lab, 
            aes(x = 0, y = 24, label = method), 
            inherit.aes = FALSE, 
            size = 2, 
            color = "black", 
            hjust = 0) 


#### Figure 3b-e ####
### Instructions 
## Take figuresofmerit files and perform heuristics and rbind into a single dataframe
# First, combine Seer NPA and NPB files
# Second, import everything minus olink and rbind
# Third, Import Olink and adjust and rbind
# Fourth, combine all quant files together for plotting

counts <- data.frame(method = character(),
                     total_peptides = integer(),
                     lod_peptides = integer(),
                     loq_peptides = integer())

## Seer NPA and NPB Combination
figures_of_merit_NPA <- fread("data/diann_output/exp4/NPA/figuresofmerit.csv") %>% mutate(NP = "NPA")
figures_of_merit_NPB <- fread("data/diann_output/exp4/NPB/figuresofmerit.csv") %>% mutate(NP = "NPB")

#add NPA and NPB column and make one object
figures_of_merit_1 <- rbind(figures_of_merit_NPA, figures_of_merit_NPB)

figures_of_merit_1 <- figures_of_merit_1 %>%
  group_by(peptide) %>%
  arrange(LOD, LOQ) %>%
  filter(row_number() == 1)

total_peptides <- nrow(figures_of_merit_1)
lod_peptides <- sum(is.infinite(figures_of_merit_1$LOD) == F)
loq_peptides <- sum(is.infinite(figures_of_merit_1$LOQ) == F)

# Manage infinite LOQ first then combine
figures_of_merit_1 <- figures_of_merit_1 %>%
  group_by(LOD, LOQ, stndev_noise, slope_linear) %>%
  mutate(LOQ = case_when(
    is.infinite(LOD) ~ Inf,
    is.infinite(LOD) == F & is.infinite(LOQ) ~ 100,
    is.infinite(LOD) == F & is.infinite(LOQ) == F ~ LOQ)) %>%
  ungroup()

figures_of_merit_1 <- figures_of_merit_1 %>%
  group_by(peptide) %>%
  arrange(LOD, LOQ) %>%
  filter(row_number() == 1)


fwrite(figures_of_merit_1, "data/diann_output/exp4/NPA_NPB_Combined/figuresofmerit.csv",
       row.names = F)

c <- data.frame(method = "seer",
                total_peptides = total_peptides,
                lod_peptides = lod_peptides,
                loq_peptides = loq_peptides)

counts <- rbind(counts, c)


## One long figures of merit dataframe 

methods <- c("neat", "acid", "preomics", "magnet")

fom_all <- data.frame(peptide = character(),
                      LOD = numeric(),
                      LOQ = numeric(),
                      slope_linear = numeric(),
                      intercept_linear = numeric(),
                      intercept_noise = numeric(),
                      stndev_noise = numeric(),
                      method = character(),
                      name = character(),
                      value = numeric())

for (i in methods) {
  
  figures_of_merit <- fread(paste0("data/diann_output/exp4/", i, "/figuresofmerit.csv")) %>% 
    mutate(method = i) 
  
  total_peptides <- nrow(figures_of_merit)
  lod_peptides <- sum(is.infinite(figures_of_merit$LOD) == F)
  loq_peptides <- sum(is.infinite(figures_of_merit$LOQ) == F)
  
  # Manage Infinite Values
  figures_of_merit <- figures_of_merit %>%
    group_by(LOD, LOQ, stndev_noise, slope_linear) %>%
    mutate(LOQ = case_when(
      is.infinite(LOD) ~ Inf,
      is.infinite(LOD) == F & is.infinite(LOQ) ~ 100,
      is.infinite(LOD) == F & is.infinite(LOQ) == F ~ LOQ)) %>%
    ungroup()
  
  # Basic EDA
  figures_of_merit_LOD <- figures_of_merit %>%
    arrange(LOD) %>%
    mutate(LOD_Percentile = percent_rank(LOD),
           LOD_Count = row_number()) %>%
    pivot_longer(cols = starts_with('LOD_'), names_to = 'name', values_to = 'value')
  
  figures_of_merit_LOQ <- figures_of_merit %>%
    arrange(LOQ) %>%
    mutate(LOQ_Percentile = percent_rank(LOQ),
           LOQ_Count = row_number()) %>%
    pivot_longer(cols = starts_with('LOQ_'), names_to = 'name', values_to = 'value')
  
  fom_all <- rbind(fom_all, figures_of_merit_LOD)
  fom_all <- rbind(fom_all, figures_of_merit_LOQ)
  
  c <- data.frame(method = i,
                  total_peptides = total_peptides,
                  lod_peptides = lod_peptides,
                  loq_peptides = loq_peptides)
  
  counts <- rbind(counts, c)
  
}

# Add Seer
figures_of_merit_seer <- fread("data/diann_output/exp4/NPA_NPB_Combined/figuresofmerit.csv",
                               data.table = F) %>%
  select(-NP) %>%
  mutate(method = "seer")

figures_of_merit_LOD <- figures_of_merit_seer %>%
  arrange(LOD) %>%
  mutate(LOD_Percentile = percent_rank(LOD),
         LOD_Count = row_number()) %>%
  pivot_longer(cols = starts_with('LOD_'), names_to = 'name', values_to = 'value')

figures_of_merit_LOQ <- figures_of_merit_seer %>%
  arrange(LOQ) %>%
  mutate(LOQ_Percentile = percent_rank(LOQ),
         LOQ_Count = row_number()) %>%
  pivot_longer(cols = starts_with('LOQ_'), names_to = 'name', values_to = 'value')

fom_all <- rbind(fom_all, figures_of_merit_LOD)
fom_all <- rbind(fom_all, figures_of_merit_LOQ)


## Add Olink too
olink <- read_parquet('data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet') %>% #Olink data
  filter(grepl('O4', SampleID) | grepl('CONTROL', SampleType))
lod <- olink_lod(olink)
olink <- olink %>% 
  left_join(lod) %>%
  filter(grepl('O4', SampleID) | grepl('NEGATIVE_CONTROL', SampleType)) %>%
  mutate(olink_suggested_signal = PCNormalizedNPX > PCNormalizedLOD) %>%
  mutate(Group = substr(SampleID, 4,4)) %>%
  mutate(concentration = case_when(
    Group == 'a' ~ 100,
    Group == 'b' ~ 70,
    Group == 'c' ~ 50,
    Group == 'd' ~ 30,
    Group == 'e' ~ 10,
    Group == 'f' ~ 5,
    Group == 'g' ~ 1,
    Group == 'h' ~ .5,
    Group == 'i' ~ .1,
    Group == 'j' ~ 0
  )) 

figures_of_merit <- fread("data/diann_output/exp4/human/figuresofmerit.csv",
                          data.table = F)

figures_of_merit <- olink %>%
  filter(SampleType != 'NEGATIVE_CONTROL') %>%
  filter(!grepl('Extension', Assay)) %>%
  group_by(Assay) %>%
  reframe(cor = cor(PCNormalizedNPX, concentration)) %>%
  mutate(peptide = Assay) %>%
  right_join(figures_of_merit)

total_peptides <- sum(!is.na(figures_of_merit$LOD))

figures_of_merit <- figures_of_merit %>%
  group_by(LOD, LOQ, stndev_noise, slope_linear, cor) %>%
  mutate(LOD = if_else(cor < 0.6, Inf, LOD),
         LOQ = if_else(cor < 0.6, Inf, LOQ)) %>%
  mutate(LOQ = case_when(
    is.infinite(LOD) ~ Inf,
    is.infinite(LOD) == F & is.infinite(LOQ) ~ 100,
    is.infinite(LOD) == F & is.infinite(LOQ) == F ~ LOQ)) %>%
  ungroup()

lod_peptides <- sum(is.infinite(figures_of_merit$LOD) == F, na.rm = T)
loq_peptides <- sum(figures_of_merit$LOQ != Inf & figures_of_merit$LOQ != 100, na.rm = T)

figures_of_merit <- figures_of_merit %>%
  select(-Assay, -cor) %>%
  mutate(method = "olink")

figures_of_merit_LOD <- figures_of_merit %>%
  arrange(LOD) %>%
  mutate(LOD_Percentile = percent_rank(LOD),
         LOD_Count = row_number()) %>%
  pivot_longer(cols = starts_with('LOD_'), names_to = 'name', values_to = 'value')

figures_of_merit_LOQ <- figures_of_merit %>%
  arrange(LOQ) %>%
  mutate(LOQ_Percentile = percent_rank(LOQ),
         LOQ_Count = row_number()) %>%
  pivot_longer(cols = starts_with('LOQ_'), names_to = 'name', values_to = 'value')

fom_all <- rbind(fom_all, figures_of_merit_LOD)
fom_all <- rbind(fom_all, figures_of_merit_LOQ)

c <- data.frame(method = "olink",
                total_peptides = total_peptides,
                lod_peptides = lod_peptides,
                loq_peptides = loq_peptides)

counts <- rbind(counts, c)


fom_all$method <- factor(fom_all$method, 
                         levels = c("neat", "acid", "preomics", "magnet", "seer", "olink"))

fwrite(fom_all, "data/processed/Exp4_FiguresOfMerit_LODLOQ_Combined_AllMethods_AfterFiltering.csv")



# calculate number of assays in each range for each method
summaLOD <- fom_all %>%
  filter(name == "LOD_Count") %>%
  group_by(method) %>%
  mutate(
    range = case_when(
      LOD < 5 ~ "<5",
      LOD >= 5 & LOD < 20 ~ "5-20",
      LOD >= 20 & LOD < 50 ~ "20-50",
      LOD >= 50 & LOD <= 99.9 ~ "50-99.9",
      LOD == Inf ~ "Inf"
    )
  ) %>%
  group_by(range, method) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(.)) %>% # remove Olink NA Values from control assays
  filter(range != Inf) # remove inf from plot


summaLOD$range <- factor(summaLOD$range, levels = c("50-99.9", "20-50", "5-20", "<5"))

p1 <- ggplot(summaLOD, aes(x = method, y = count, fill = method, alpha = range)) + 
  geom_bar(stat = "identity", 
           position = "stack", 
           color = NA,
           width = 0.4) + 
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            hjust = 1,
            size = 2, 
            color = "black") +
  scale_fill_manual(values = pal) +
  scale_alpha_manual(values = c(0.4, 0.6, 0.8, 1)) +
  xlab(NULL) +
  ylab("# of Proteins") +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = "none") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "cm")
  ) 


summaLOQ <- fom_all %>%
  filter(name == "LOQ_Count") %>%
  group_by(method) %>%
  mutate(
    range = case_when(
      LOQ < 5 ~ "<5",
      LOQ >= 5 & LOQ < 20 ~ "5-20",
      LOQ >= 20 & LOQ < 50 ~ "20-50",
      LOQ >= 50 & LOQ <= 99.9 ~ "50-99.9",
      LOQ == 100 ~ "100",
      LOD == Inf ~ "Inf"
    )
  ) %>%
  group_by(range, method) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(.)) %>% # remove Olink NA Values from control assays
  filter(range != Inf) # remove inf from plot

summaLOQ$range <- factor(summaLOQ$range, levels = c("100", "50-99.9", "20-50", "5-20", "<5"))

p2 <- ggplot(summaLOQ, aes(x = method, y = count, fill = method, alpha = range)) + 
  geom_bar(stat = "identity", 
           position = "stack", 
           color = NA,
           width = 0.4) + 
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            hjust = 1,
            size = 2, 
            color = "black") +
  scale_fill_manual(values = pal) +
  scale_alpha_manual(values = c(0.2, 0.4, 0.6, 0.8, 1)) +
  xlab(NULL) +
  ylab("# of Proteins") +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = "none") +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "cm")
  ) 

p1 + p2

ggsave("reports/figures/Exp4_ProportionMeasuredinLODQbins_small_changedbins.pdf", 
       width = 16, height = 6, units = "cm")



## Detected or Quantified LOD or LOQ

counts1 <- counts %>%
  mutate(lod_noquantify = total_peptides - lod_peptides) %>%
  mutate(loq_noquantify = total_peptides - loq_peptides)

counts_long <- counts1 %>%
  pivot_longer(cols = -c(total_peptides, method), 
               names_to = "Type", 
               values_to = "Value") %>%
  select(-total_peptides) %>%
  separate(Type, c("lodq", "type")) %>%
  mutate(method = factor(method, levels = c("neat", "acid", "preomics", "magnet", "seer", "olink")))


p3 <- ggplot(counts_long %>% filter(lodq == "lod"), aes(x = method, y = Value, alpha = type)) + 
  geom_bar(stat = "identity", 
           position = "stack", 
           color = NA,
           width = 0.4,
           fill = "darkgray") + 
  geom_text(aes(label = Value), 
            position = position_stack(vjust = 0.5), 
            hjust = 1,
            size = 2, 
            color = "black") +
  scale_alpha_manual(values = c(0.6, 1)) +
  xlab(NULL) +
  ylab("# of Proteins") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "cm")
  ) 

p4 <- ggplot(counts_long %>% filter(lodq == "loq"), aes(x = method, y = Value, alpha = type)) + 
  geom_bar(stat = "identity", 
           position = "stack", 
           color = NA,
           width = 0.4,
           fill = "darkgray") + 
  geom_text(aes(label = Value), 
            position = position_stack(vjust = 0.5), 
            hjust = 1,
            size = 2, 
            color = "black") +
  scale_alpha_manual(values = c(0.6, 1)) +
  xlab(NULL) +
  ylab("# of Proteins") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.2, "cm")
  )  

p3 + p4