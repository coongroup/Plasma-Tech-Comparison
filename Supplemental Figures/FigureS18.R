# Copyright 2025 Seer, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

################################################################################

### OLink Data Analysis
rm(list = ls())
library(tidyverse)
library(data.table)
library(OlinkAnalyze)
library(arrow)

olink <- read_parquet('data/olink_output/P2024-198-OLI_HT_93_NPX_2024-10-10.parquet') %>% #Olink data
  filter(grepl('O4', SampleID) | SampleType == 'NEGATIVE_CONTROL')

lod <- olink_lod(olink)

olink <- olink %>% 
  left_join(lod) %>%
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

figures_of_merit <- fread('data/olink_output/human/figuresofmerit.csv') #human figures of merit file
figures_of_merit <- olink %>%
  filter(SampleType != 'NEGATIVE_CONTROL') %>%
  filter(!grepl('Extension', Assay)) %>%
  group_by(Assay) %>%
  reframe(cor = cor(PCNormalizedNPX, concentration)) %>%
  dplyr::rename(peptide = Assay) %>%
  right_join(figures_of_merit)
figures_of_merit <- figures_of_merit %>%
  group_by(LOD, LOQ, stndev_noise, slope_linear, cor) %>%
  mutate(LOD = case_when(
    is.infinite(stndev_noise) & slope_linear <= 0.000 & is.infinite(LOD) ~ Inf, # not quantitative
    is.infinite(stndev_noise) & slope_linear > 0 & is.infinite(LOD) & cor >= 0.7 ~ 0, # infinitely quantitative to assay
    is.infinite(stndev_noise) & slope_linear > 0 & is.infinite(LOD) & cor < 0.7 ~ Inf, # infinitely non-quantitative to assay, olink mod
    is.infinite(LOD) == F ~ LOD # quantitative to assay
  )) %>%
  mutate(LOQ = case_when(
    is.infinite(LOD) ~ Inf,
    is.infinite(LOD) == F & is.infinite(LOQ) & slope_linear > 0.001 ~ LOD,
    is.infinite(LOD) == F & is.infinite(LOQ) == F ~ LOQ,
    LOD == 0 & is.infinite(LOQ) & slope_linear > 0.000 ~ 0,
    LOD == 0 & is.infinite(LOQ) & slope_linear <= 0 ~ Inf)) %>%
  ungroup()

figures_of_merit2 <- fread('data/olink_output/chicken/figuresofmerit.csv') #chicken figures of merit file
figures_of_merit2 <- olink %>%
  filter(SampleType != 'NEGATIVE_CONTROL') %>%
  filter(!grepl('Extension', Assay)) %>%
  group_by(Assay) %>%
  reframe(cor = cor(PCNormalizedNPX, 100-concentration)) %>%
  dplyr::rename(peptide = Assay) %>%
  right_join(figures_of_merit2)
figures_of_merit2 <- figures_of_merit2 %>%
  group_by(LOD, LOQ, stndev_noise, slope_linear, cor) %>%
  mutate(LOD = case_when(
    is.infinite(stndev_noise) & slope_linear <= 0.000 & is.infinite(LOD) ~ Inf, # not quantitative
    is.infinite(stndev_noise) & slope_linear > 0 & is.infinite(LOD) & cor >= 0.6 ~ 0, # infinitely quantitative to assay
    is.infinite(stndev_noise) & slope_linear > 0 & is.infinite(LOD) & cor < 0.6 ~ Inf, # infinitely non-quantitative to assay, olink mod
    is.infinite(LOD) == F ~ LOD # quantitative to assay
  )) %>%
  mutate(LOQ = case_when(
    is.infinite(LOD) ~ Inf,
    is.infinite(LOD) == F & is.infinite(LOQ) & slope_linear > 0.001 ~ LOD,
    is.infinite(LOD) == F & is.infinite(LOQ) == F ~ LOQ,
    LOD == 0 & is.infinite(LOQ) & slope_linear > 0.000 ~ 0,
    LOD == 0 & is.infinite(LOQ) & slope_linear <= 0 ~ Inf)) %>%
  ungroup()

human_lod_peptides <- figures_of_merit %>% filter(!is.infinite(LOD)) %>% pull(peptide)
chicken_lod_peptides <- figures_of_merit2 %>% filter(!is.infinite(LOD)) %>% pull(peptide)

#### Figure S18 ####
olink %>%
  filter(!grepl('control', Assay)) %>%
  group_by(Assay) %>%
  mutate(mean_neg = mean(2**PCNormalizedNPX[SampleType == 'NEGATIVE_CONTROL'], na.rm = T),
         sd_neg = sd(2**PCNormalizedNPX[SampleType == 'NEGATIVE_CONTROL'], na.rm = T),
         mad_neg = mad(2**PCNormalizedNPX[SampleType == 'NEGATIVE_CONTROL'], na.rm = T),
         measurable = 2**PCNormalizedNPX > mean_neg + 3*sd_neg,
         measurable_robust = 2**PCNormalizedNPX > mean_neg + 3*mad_neg) %>%
  filter(!grepl('CONTROL', SampleType)) %>%
  pivot_longer(cols = c(olink_suggested_signal, measurable, measurable_robust)) %>%
  group_by(Assay) %>%
  mutate(MMCC_LOD_Assigned = case_when(
    Assay %in% human_lod_peptides ~ 'HUMAN Response',
    Assay %in% chicken_lod_peptides ~ 'CHICK Response',
    Assay %in% chicken_lod_peptides & Assay %in% human_lod_peptides ~ 'Not Measurable',
    Assay %in% chicken_lod_peptides == F & Assay %in% human_lod_peptides == F ~ 'Not Measurable'
  )) %>%
  group_by(concentration, Assay, name, MMCC_LOD_Assigned) %>%
  reframe(proportion_measured = sum(value) / n()) %>%
  mutate(proportion_measured = (round(3*proportion_measured)/3)) %>%
  mutate(prop_measured_label = round(proportion_measured, 2)) %>%
  filter(name == 'measurable') %>% #Toggle to make data easier to interpret 
  mutate(name = case_when(
    name == 'measurable' ~ 'Mean + 3*StDev.',
    name == 'measurable_robust' ~ 'Mean + 3*MAD',
    name == 'olink_suggested_signal' ~ 'OLink Analyze Method'
  )) %>%
  ggplot(aes(factor(prop_measured_label), fill = factor(concentration))) + 
  geom_bar(position = 'dodge', col = 'black') + 
  theme_light(base_size = 16) + facet_grid(MMCC_LOD_Assigned~name, scales = 'free_y') + 
  labs(x = 'Measurement Reproducibility', 
       y = 'Assays Above Threshold at Concentration')  + theme(legend.position = 'bottom') + 
  guides(fill = guide_legend(title = 'Assay Concentration %', ncol = 6))

