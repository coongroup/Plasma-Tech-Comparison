
library(devtools)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(ggpmisc)



#read in all 6 reports
file_a <- "20240918_WFB_PlasmaTechComparison_Experiment1_acidt_dDIA_Report_WFB_Report (Normal).tsv"
file_n <- "20240918_WFB_PlasmaTechComparison_Experiment1_neat_dDIA_Report_WFB_Report (Normal).tsv"
file_p <- "20240918_WFB_PlasmaTechComparison_Experiment1_preomics_dDIA_Report_WFB_Report (Normal).tsv"
file_m <- "20240918_WFB_PlasmaTechComparison_Experiment1_magnet_dDIA_Report_WFB_Report (Normal).tsv"
file_s <- "20240918_WFB_PlasmaTechComparison_Experiment1_seer_dDIA_Report_WFB_Report (Normal).tsv"

files <- list(neat = file_n, 
              acid = file_a, 
              preomics = file_p, 
              magnet = file_m, 
              seer = file_s)

for (i in names(files)) {
  #read in file
  spec <- read.delim(paste0("data/spectronaut_output/Experiment1_separate/", files[[i]]), sep = "\t")
  
  #Pivot report into protein quant values for each sample (can adjust parameters to do peptide too)
  selects <- (is.na(spec$PG.Qvalue) | spec$PG.Qvalue <= 0.01) &
    (is.na(spec$EG.Qvalue) | spec$EG.Qvalue <= 0.01) 
  spec_raw <- spec[selects, c("R.FileName", "PEP.NrOfMissedCleavages", "EG.ModifiedSequence", "EG.TotalQuantity..Settings.")]
  spec_raw %>%
    dplyr::group_by(EG.ModifiedSequence, PEP.NrOfMissedCleavages) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n > 1L)  
  spec_pep <- pivot_wider(spec_raw,
                          names_from = R.FileName,
                          values_from = EG.TotalQuantity..Settings.,
                          values_fill = NA,
                          values_fn = list(EG.TotalQuantity..Settings. = mean))
  spec_pep <- spec_pep %>%
    group_by(EG.ModifiedSequence) %>%
    summarise_all(~ na.omit(.) %>% .[1]) %>%
    ungroup()
  
  #make separate dataframe
  df_name <- paste0(i, "_pep")
  assign(df_name, spec_pep)
  
}


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



#### Figure S5a ####
dfs <- list(neat_pep = neat_pep, 
            acid_pep = acid_pep, 
            magnet_pep = magnet_pep, 
            preomics_pep = preomics_pep, 
            seer_pep = seer_pep)

digest_all <- data.frame(NumMissCleav = c(0, 1, 2))

for (i in names(dfs)) {
  
  y <- dfs[[i]]
  digest <- data.frame(table(y$PEP.NrOfMissedCleavages))
  
  
  colnames(digest) <- c("NumMissCleav", paste0(i, "_count"))
  
  # make combined
  digest_all <- merge(digest_all, digest, by = "NumMissCleav")
  
}

colnames(digest_all) <- c("NumMissCleav", "Neat", "Acid", "MagNet", "Preomics", "Seer")

digest_long <- digest_all %>%
  pivot_longer(cols = !contains("NumMissCleav"), 
               names_to = "Method", 
               values_to = "Value")
digest_long <- digest_long %>%
  group_by(Method) %>%
  mutate(Proportion = Value / sum(Value) * 100)
digest_long <- digest_long %>%
  group_by(Method) %>%
  mutate(NumMissCleav = factor(NumMissCleav, levels = c(2, 1, 0)))
digest_long$Method <- factor(digest_long$Method, levels = c("Neat", "Acid", "Preomics", "MagNet", "Seer"))


#make plot
ggplot(digest_long, aes(x = Method, y = Value, fill = NumMissCleav)) + 
  geom_bar(stat = "identity", 
           position = "stack", 
           color = NA) + 
  geom_text(aes(label = paste0(round(Proportion, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 2, 
            color = "black") +
  scale_fill_manual(values = brewer.pal(3, "YlOrBr")) +
  xlab(NULL) +
  ylab("# of Peptides") +
  scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(expand = c(0,0), breaks = seq(0, 40, by = 8)) +
  guides(fill = guide_legend(title = "NumMissCleav")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.2, "in")
  ) 



#### Figure S5b ####
#Look at proportion of unspecific/semi-specific for each method
dfs <- list(neat_pep = neat_pep, acid_pep = acid_pep, magnet_pep = magnet_pep, preomics_pep = preomics_pep, seer_pep = seer_pep)

digest_all <- data.frame(type = c("Specific", "Specific-C", "Specific-N", "Unspecific"))

for (i in names(dfs)) {
  
  y <- dfs[[i]]
  digest <- data.frame(table(y$PEP.DigestType....Trypsin.P.))
  
  digest <- digest %>%
    filter(!grepl(",", Var1))
  
  colnames(digest) <- c("type", paste0(i, "_count"))
  
  # make combined
  digest_all <- merge(digest_all, digest, by = "type")
  
}

colnames(digest_all) <- c("type", "Neat", "Acid", "MagNet", "Preomics", "Seer")

digest_long <- digest_all %>%
  pivot_longer(cols = !contains("type"), 
               names_to = "Method", 
               values_to = "Value")
digest_long <- digest_long %>%
  group_by(Method) %>%
  mutate(Proportion = Value / sum(Value) * 100)
digest_long <- digest_long %>%
  group_by(Method) %>%
  mutate(type = factor(type, levels = c("Unspecific", "Specific-C", "Specific-N", "Specific")))
digest_long$Method <- factor(digest_long$Method, levels = c("Neat", "Acid", "Preomics", "MagNet", "Seer"))

#make plot
ggplot(digest_long, aes(x = Method, y = Value, fill = type)) + 
  geom_bar(stat = "identity", position = "stack", color = NA) + 
  geom_text(aes(label = paste0(round(Proportion, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 2, color = "black") +
  scale_fill_manual(values = brewer.pal(4, "YlOrBr")) +
  #ggtitle("Unspecific Digest") +
  xlab(NULL) +
  ylab("# of Peptides") +
  #ylim(-10, 0) +
  scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(expand = c(0,0), breaks = seq(0, 40, by = 8)) +
  guides(fill = guide_legend(title="Type")) +
  theme_classic() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"),
        legend.margin = margin(2, 2, 2, 2),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.2, "in")
  ) 

