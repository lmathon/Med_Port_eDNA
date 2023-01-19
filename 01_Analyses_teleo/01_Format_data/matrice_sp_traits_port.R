## Librairies
library(dplyr)
library(rsq)
library(fastDummies)
library(tidyr)
library(tidyverse)
library(purrr)
library(scales)
library(margins)
library(ggplot2)
library(gridExtra)
library(viridis)
library(ggpubr)
library(fishtree)


## Load data
data <- read.csv("01_Analyses_teleo/00_data/matrice_teleo_port.csv", header=T, row.names=1)
meta <- read.csv("00_Metadata/metadata_port.csv", header=T)
traits <- read.csv("01_Analyses_teleo/00_data/Functional_data_corrected_20220124.csv", header=T)
data$Species <- sub(pattern = "_",  replacement = " ", data$Species)
data2 <- data 
data <- column_to_rownames(data, var = "Species")

# list of species
species <- rownames(data)
# list of samples
samples <- colnames(data)

#########################################################################
## trouver les espèces en danger et exotiques
data_combined <- as.data.frame(t(data)) %>%
  rownames_to_column(var="code_spygen") %>%
  left_join(meta, by="code_spygen")

traits <- traits %>%
  filter(Species %in% species)

## Endangered species
iucn_sp <- traits %>%
  filter(IUCN_Red_List_Category %in% c("CR", "EN", "VU")) %>%
  select(Species, IUCN_Red_List_Category)

data_combined %>%
  select(any_of(c("code_spygen", "station", iucn_sp$Species)))

## Exotic species
exo_sp <- traits %>%
  filter(NI == 1) %>%
  select(Species, NI)

data_combined %>%
  select(any_of(c("code_spygen", "station", exo_sp$Species)))

## elasmobranches
elasm_sp <- traits %>%
  filter(SHarK == 1) %>%
  select(Species, SHarK)

data_combined %>%
  select(any_of(c("code_spygen", "station", elasm_sp$Species)))

## commercial
commercial_sp <- traits %>%
  filter(all_commercial_level == 1) %>%
  select(Species, all_commercial_level)

data_combined %>%
  select(any_of(c("code_spygen", "station", commercial_sp$Species)))

## pelagic
pelagic_sp <- traits %>%
  filter(Vertical_Distribution == "Pelagic") %>%
  select(Species, Vertical_Distribution)

data_combined %>%
  select(any_of(c("code_spygen", "station", pelagic_sp$Species)))


## matrix results 
sp_port_tot_traits <- full_join(data2["Species"],commercial_sp, exo_sp, by = "Species")
sp_port_tot_traits <- full_join(sp_port_tot_traits, exo_sp, by = "Species")
sp_port_tot_traits <- full_join(sp_port_tot_traits, iucn_sp, by = "Species") 
sp_port_tot_traits <- full_join(sp_port_tot_traits, pelagic_sp, by = "Species") %>%
  replace(is.na(.), 0)

write.csv(sp_port_tot_traits, "01_Analyses_teleo/00_data/sp_port_tot_traits.csv")

## liste des espèces qui sont pas dans la matrice de traits du tout 
traits <- read.csv("01_Analyses_teleo/00_data/Functional_data_corrected_20220124.csv", header=T) 

sp_port_not_in_trait <- full_join(data2, traits, by = "Species")
sp_port_not_in_trait <- sp_port_not_in_trait[,1:19]

sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]
sp_port_not_in_trait <- sp_port_not_in_trait[,-2]


