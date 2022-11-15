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
data <- read.csv2("teleo_presence.csv")

meta <- read.csv("metadata_eDNA_port_october_2021.csv", header=T)

traits <- read.csv("C:/Users/AliciaDalongeville/Documents/MARBEC/LABCOM/Data/ADNe/Functional_data_corrected_20220124.csv", header=T)


##########################################################################################
###### Calculate indicators 
##########################################################################################
library(dplyr)
library(tidyr)
library(fastDummies)
library(ape)
library(fishtree)
library(picante)
library(geiger)
loadNamespace("rfishbase")

# make species names row names
data <- data %>%
  column_to_rownames(var="scientific_name")

# list of species
species <- rownames(data)
# list of samples
samples <- colnames(data)

# Quelles espèces ne sont pas dans les traits ?
#check <- species[which(species %in% traits$Species ==F)]
#print(check)
# Remove the 12 species from the data for now : il faudra les checker manuellement
#data2 <- data %>%
#  filter(rownames(.) %in% check ==F)
#write.csv(data2, file="SC21250-CEFE_teleo_presence_refined.csv")

#######################################################################################################
## Calculate the indicators for each sample
#######################################################################################################
## Create the result matrix
indicators <- matrix(NA,ncol(data), 14,
                     dimnames=list(colnames(data),
                                   c("R",  "FD", "LFI", "Crypto", "DeBRa",  "CTI", "Exo","RedList", "Chondri", "Trophic", "Commercial", "High_commerc", "PD", "Vulner")))

##########################################################
## 1 - Species Richness R
##########################################################
indicators[,1] <- apply(data, 2, sum)

##########################################################
## 2 - Functional diversity FD
##########################################################
for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- data %>%
    dplyr::filter(data[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Get the functional groups of these species
  fd_i <- as.factor(traits[which(traits$Species %in% s_i), "m.cluster_core"])
  
  # Number of unique functional groups
  indicators[i,2] <- length(unique(fd_i))
}

##########################################################
## 3 - Large Reef Fish Index - LFI
## 5 - Ratio Demerso-pelagic / benthic
## 7 - Non-indigenous species
## 9 - Chondrichtyen species
## 11 - Commercial species
## 12 - Highly commercial species 
##########################################################
for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- data %>%
    filter(data[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Calculate indicators
  indicators[i,"LFI"] <- sum(traits[which(traits$Species %in% s_i), "LRFI"]) 
  indicators[i,"DeBRa"] <- sum(traits[which(traits$Species %in% s_i), "DP"]) / (sum(traits[which(traits$Species %in% s_i), "B"])+1)
  indicators[i,"Exo"] <- sum(traits[which(traits$Species %in% s_i), "NI"])
  indicators[i,"Chondri"] <- sum(traits[which(traits$Species %in% s_i), "SHarK"])
  indicators[i,"Commercial"] <- sum(traits[which(traits$Species %in% s_i), "all_commercial_level"])
  indicators[i,"High_commerc"] <- sum(traits[which(traits$Species %in% s_i), "highly_commercial_only"])
}

###########################################################
## 4 Cryptobenthic (definition Brandl et al. 2018)
# Brandl SJ, Goatley CHR, Bellwood DR, Tornabene L. 2018 The hidden half: ecology and evolution of cryptobenthic fishes on coral reefs. Biol. Rev. 93, . (doi:10.1111/brv.12423))
###########################################################
crypto_families = c("Tripterygiidae", "Grammatidae", "Creediidae", "Aploactinidae", "Gobiidae", 
                    "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", 
                    "Plesiopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae")

traits <- traits %>%
  mutate(crypto_Brandl = if_else(Family %in% crypto_families, 1,0))

for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- data %>%
    filter(data[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Calculate indicator
  indicators[i,"Crypto"] <- sum(traits[which(traits$Species %in% s_i), "crypto_Brandl"]) 
}

##########################################################
## 6 - Community Temperature Index - CTI
## 10 - Mean trophic level
## 14 - Vulnerability
##########################################################
for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- data %>%
    filter(data[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # calculate indicators
  indicators[i,"CTI"] <- mean(traits[which(traits$Species %in% s_i), "SSTmean"], na.rm = T) 
  indicators[i,"Trophic"] <- mean(traits[which(traits$Species %in% s_i), "Trophic_level"], na.rm=T) 
  indicators[i,"Vulner"] <- mean(traits[which(traits$Species %in% s_i), "Vulnerability"], na.rm=T)
}

##########################################################
## 8 - Red List IUCN
##########################################################
## Assign different weighs to the IUCN categories: NT/LC = 0,VU = 1, EN = 2, CR = 3

### Make a dummy variable for IUCN categories
traits <- dummy_cols(traits, select_columns = 'IUCN_Red_List_Category')

## Calculate indicator
for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- data %>%
    filter(data[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Number of species per category
  VU <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_VU"], na.rm=T)
  EN <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_EN"], na.rm=T)
  CR <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_CR"], na.rm=T)
  
  # Calculate indicator
  indicators[i,"RedList"] <- 0 + VU + EN*2 + CR*3
}

##########################################################
## 13 - Phylogenetic Diversity - PD
##########################################################
# Retrieve the phylogeny of only native reef species across all three oceans.
phy <- fishtree_phylogeny(species = species)

plot(phy, show.tip.label = FALSE)
tiplabels(tip = which(phy$tip.label %in% species),
          pch=19, cex=2)

rownames(data) <- gsub(" ", "_", rownames(data), fixed = TRUE)

## check that phylogeny and data have matching names
#nc <- geiger::name.check(phy, data) # 35 species not in tree
#missing <- gsub("_", " ", nc$data_not_tree, fixed = TRUE)
## List of missing species that are not Chondrichtyens (as they are not in the tree)
#missing_fish <- traits[which(traits$Species %in% missing & traits$Super_class != "chondrichtyen"),"Species"]

# Manually check synonyms and find the species
species[species == "Mullus barbatus"] <- "Mullus barbatus barbatus"
rownames(data)[rownames(data) == "Mullus_barbatus"] <- "Mullus_barbatus_barbatus"

species[species == "Diplodus sargus"] <- "Diplodus sargus sargus"
rownames(data)[rownames(data) == "Diplodus_sargus"] <- "Diplodus_sargus_sargus"

species[species == "Diplodus cervinus"] <- "Diplodus cervinus cervinus"
rownames(data)[rownames(data) == "Diplodus_cervinus"] <- "Diplodus_cervinus_cervinus"

species[species == "Chelon auratus"] <- "Liza aurata"
rownames(data)[rownames(data) == "Chelon_auratus"] <- "Liza_aurata"

species[species == "Chelon ramada"] <- "Liza ramada"
rownames(data)[rownames(data) == "Chelon_ramada"] <- "Liza_ramada"

# Retrieve the missing phylogeny 
phy <- fishtree_phylogeny(species = species, type="phylogram")
nc <- geiger::name.check(phy, data) # 34 espèces manquantes

# Remove from the data the species that are not in the tree
data2 <- data[which(rownames(data) %in% nc$data_not_tree == F),]

# Transpose the ADNe matrix 
data2 <- t(data2)

# prune the tree
prunedTree <- prune.sample(data2,phy)

# Calculate PD
pd.result <- pd(data2, prunedTree, include.root=T)

# Add PD to indicator dataframe
indicators[,"PD"] <- pd.result$PD

write.table(indicators, file="Indicators_ports_2021.csv", sep=",")

#### Summarize indicators per site
indicators  %>% 
  as.data.frame(.) %>%
  rownames_to_column(var="code_spygen") %>%
  inner_join(meta, by="code_spygen") %>%
  select(R, Crypto, Chondri, Commercial, Exo, RedList, site) %>%
  group_by(site) %>%
  summarise(across(everything(), list(mean)))

#### Summarize indicators per habitat (biohut ou non)
indicators  %>% 
  as.data.frame(.) %>%
  rownames_to_column(var="code_spygen") %>%
  inner_join(meta, by="code_spygen") %>%
  # selectionner les sites où les biohuts on été échantillonné
  filter(site %in% c("Cannes", "Marseillan ", "Saintes Maries de la Mer")) %>%
  select(R, Crypto, Chondri, Commercial, Exo, RedList, habitat) %>%
  group_by(habitat) %>%
  summarise(across(everything(), list(mean)))

#### Summarize total
indicators  %>% 
  as.data.frame(.) %>%
  summarise(across(everything(), list(mean=mean,
                                      min=min,
                                      max=max)))

#######################################################################################################
## Plots
#######################################################################################################
toplot <- indicators %>%
  as.data.frame(.) %>%
  rownames_to_column(var="code_spygen") %>%
  inner_join(meta, by="code_spygen") %>%
  # selectionner les sites où les biohuts on été échantillonné
  filter(site %in% c("Cannes", "Marseillan ", "Saintes Maries de la Mer")) %>%
  dplyr::select(code_spygen, R, Crypto, Chondri, Commercial, Exo, 
                RedList, habitat, site)
colnames(toplot)[c(2:7)] <- c("Richesse spécifique", "Cryptobenthiques", 
                                "Requins & raies", "Espèces commerciales" ,
                              "Espèces exotiques", "Espèces menacées")

## Create plots
response = names(toplot)[2:7]
expl = names(toplot)[8]

response = set_names(response)
response

expl = set_names(expl)
expl

# Create function to make a boxplot
boxplot_fun = function(x, y) {
  ggplot(toplot, aes(x = .data[[x]], y = .data[[y]]) ) +
    geom_boxplot() +
    theme_bw() +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=14),
          axis.title.x=element_blank(),
          axis.title.y =element_text(size=15) ) +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size=4)# +
    #rotate_x_text(30)
}

## Loop tidyr to make the plots for all response and all explanatory (here 1) variables
resp_expl = tidyr::expand_grid(response, expl) # dataframe response/expla var names
resp_expl
# Make the plots
all_plots = pmap(resp_expl, ~ boxplot_fun(x = .y, y = .x) )

## Print the plots in a grid
png("boxplots_indicateurs_habitat.png", width = 900, height = 900)
cowplot::plot_grid(plotlist = all_plots)
dev.off()

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

## chondri
elasm_sp <- traits %>%
  filter(SHarK == 1) %>%
  select(Species, SHarK)

data_combined %>%
  select(any_of(c("code_spygen", "station", elasm_sp$Species)))

