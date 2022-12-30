### Load the libraries
library(dplyr)
library(rsq)
library(fastDummies)
library(tidyr)
library(purrr)
library(scales)
library(margins)
library(ggplot2)
library(viridis)
library(spdep)
library(ggpubr)
library(gridExtra)
library(spaMM)
library(DHARMa)
library(visreg)
library(multcompView)
library(multcomp)

# load data ports
ind_ports <- read.csv("01_Analyses_teleo/00_data/indicators_ports_2022_per_filter.csv", header=T, row.names=1)

# Charger datas milieu naturel (= réserve et hors réserve)
meta_nat <- read.csv2("00_Metadata/metadata_milieu_naturel.csv", header=T)
ind_nat <- read.csv("01_Analyses_teleo/00_data/indicators_milieu_naturel.csv", header=T, row.names = 1) %>%
  rownames_to_column(var="Sample") %>%
  inner_join(meta_nat, by=c("Sample"="SPYGEN_code")) %>%
  rename(habitat = protection,
         X = longitude_start_DD,
         Y = latitude_start_DD) %>%
  column_to_rownames(var="Sample") %>%
  dplyr::select(colnames(ind_ports), "habitat", "Confinement", "X", "Y") %>%
  mutate(Confinement = case_when(Confinement == 'Y' ~ 'Lockdown',
                         TRUE ~ 'Unlock'))

### dataframe ports
meta_ports <- read.csv("00_Metadata/metadata_port.csv", header=T)
ind_ports <- ind_ports %>%
  rownames_to_column(var="Sample") %>%
  inner_join(meta_ports, by=c("Sample"="code_spygen")) %>%
  mutate(Confinement = "Unlock") %>%
  column_to_rownames(var="Sample") %>%
  dplyr::select(colnames(ind_nat))

# Combiner les deux datasets
ind_all <- rbind(ind_nat,ind_ports) %>%
  mutate(DeBRa = log10(DeBRa)) %>%
  mutate_at('Confinement', as.factor)  %>%
  mutate_at('X', as.numeric)  %>%
  mutate_at('Y', as.numeric)


ind_all$habitat <- factor(ind_all$habitat, levels=c("reserve", "outside", "Port"))
ind_all$Confinement <- factor(ind_all$Confinement, labels=c("Unlock", "Lockdown"))

ind_names <- c("Total species richness", "Cryptobenthic species richness", "Threatened species richness", "Commercial species richness")  


### Matrix of spatial coordinates
coords<- ind_all %>%
  dplyr::select(X, Y) %>%
  as.matrix(.)



###################################################################
## Fit models
###################################################################
# Fit the model
mod <- list()
mod[[1]] <- fitme(R ~ habitat * Confinement + Matern(1|X+Y),
                  family = "gaussian",  data = ind_all)

mod[[2]] <- fitme(Crypto ~ habitat * Confinement + Matern(1|X+Y),
                  family = "gaussian",  data = ind_all)

mod[[3]] <- fitme(RedList ~ habitat * Confinement + Matern(1|X+Y),
                  family = "poisson",  data = ind_all)

mod[[4]] <- fitme(Commercial ~ habitat * Confinement + Matern(1|X+Y),
                  family = "gaussian",  data = ind_all)


names(mod) <- ind_names
saveRDS(mod, file="01_Analyses_teleo/02_Analyses/models_spaMM_indicators.RDS")
