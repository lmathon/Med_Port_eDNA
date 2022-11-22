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
data <- read.csv("01_Analyses_teleo/00_data/teleo_presence.csv")
data$scientific_name <- gsub("_", " ", data$scientific_name)

meta <- read.csv("00_Metadata/metadata_port.csv", header=T)
# Add a column "campaign" : October21 or June22
#meta <- meta %>%
#  mutate(Campaign = ifelse(grepl("2022", meta$date,fixed = TRUE), 'June22', 'October21'))
# save
#write.csv(meta, "00_Metadata/metadata_port.csv", row.names=F)

traits <- read.csv("01_Analyses_teleo/00_data/Functional_data_corrected_20220124.csv", header=T)

taxo <- read.csv("01_Analyses_teleo/00_data/teleo_taxo.csv")

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
# list elasmobranch species
elasmo <- taxo %>%
  filter(Class == "Chondrichthyes") %>%
  pull(scientific_name)

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
  indicators[i,"LFI"] <- sum(traits[which(traits$Species %in% s_i), "LRFI"], na.rm=T) 
  indicators[i,"DeBRa"] <- sum(traits[which(traits$Species %in% s_i), "DP"], na.rm=T) / (sum(traits[which(traits$Species %in% s_i), "B"], na.rm=T)+1)
  indicators[i,"Exo"] <- sum(traits[which(traits$Species %in% s_i), "NI"], na.rm=T)
  indicators[i,"Commercial"] <- sum(traits[which(traits$Species %in% s_i), "all_commercial_level"], na.rm=T)
  indicators[i,"High_commerc"] <- sum(traits[which(traits$Species %in% s_i), "highly_commercial_only"], na.rm=T)
}

###########################################################
## 4 Cryptobenthic (definition Brandl et al. 2018)
# Brandl SJ, Goatley CHR, Bellwood DR, Tornabene L. 2018 The hidden half: ecology and evolution of cryptobenthic fishes on coral reefs. Biol. Rev. 93, . (doi:10.1111/brv.12423))
###########################################################
crypto_families = c("Tripterygiidae", "Grammatidae", "Creediidae", "Aploactinidae", "Gobiidae", 
                    "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", 
                    "Plesiopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae")

crypto <- taxo %>%
  filter(Family %in% crypto_families) %>%
  pull(scientific_name)


for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- data %>%
    filter(data[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  crypto_i <- s_i[which(s_i %in% crypto)]
  
  # Calculate indicator
  indicators[i,"Crypto"] <- length(crypto_i)
}

###########################################################
## 9 Chondrichtyens
###########################################################
for (i in 1:nrow(indicators)) { # for each sample
  # list species present in the sample
  s_i <- data %>%
    filter(data[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  elasmo_i <- s_i[which(s_i %in% elasmo)]
  
  # Calculate indicator
  indicators[i,"Chondri"] <-length(elasmo_i)
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
  indicators[i,"RedList"] <- VU + EN + CR
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
nc <- geiger::name.check(phy, data) # 40 espèces manquantes

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

write.table(indicators, file="01_Analyses_teleo/00_data/Indicators_ports_2022_per_filter.csv", sep=",")


#### Summarize total
indicators  %>% 
  as.data.frame(.) %>%
  summarise(across(everything(), list(mean=mean,
                                      min=min,
                                      max=max)))

#######################################################################################################
## Calculate indicators per site (pooling the species list of the two replicates)
#######################################################################################################
rownames(data) <- gsub("_", " ", rownames(data))

data2 <- data %>%
  t(.) %>%
  as.data.frame(.) %>%
  rownames_to_column(var="code_spygen") %>%
  left_join(meta, by="code_spygen") %>%
  # remove the biohut samples
  filter(habitat == "Port") %>%
  group_by(site, Campaign) %>%
  summarise_at(2:132, sum) %>% 
  # convert to presence/abscence
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  # create a new var combining site and campaign column
  mutate(Names = paste(site,Campaign, sep="_")) %>%
  column_to_rownames(var="Names") %>%
  select(3:ncol(.)) %>%
  t() %>%
  as.data.frame()

## Create the result matrix
indicators2 <- matrix(NA,ncol(data2), 14,
                     dimnames=list(colnames(data2),
                                   c("R",  "FD", "LFI", "Crypto", "DeBRa",  "CTI", "Exo","RedList", "Chondri", "Trophic", "Commercial", "High_commerc", "PD", "Vulner")))

##########################################################
## 1 - Species Richness R
##########################################################
indicators2[,1] <- apply(data2, 2, sum)

##########################################################
## 2 - Functional diversity FD
##########################################################
for (i in 1:nrow(indicators2)) { # for each sample
  # list species present in the sample
  s_i <- data2 %>%
    dplyr::filter(data2[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Get the functional groups of these species
  fd_i <- as.factor(traits[which(traits$Species %in% s_i), "m.cluster_core"])
  
  # Number of unique functional groups
  indicators2[i,2] <- length(unique(fd_i))
}

##########################################################
## 3 - Large Reef Fish Index - LFI
## 5 - Ratio Demerso-pelagic / benthic
## 7 - Non-indigenous species
## 11 - Commercial species
## 12 - Highly commercial species 
##########################################################
for (i in 1:nrow(indicators2)) { # for each sample
  # list species present in the sample
  s_i <- data2 %>%
    filter(data2[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Calculate indicators2
  indicators2[i,"LFI"] <- sum(traits[which(traits$Species %in% s_i), "LRFI"], na.rm=T) 
  indicators2[i,"DeBRa"] <- sum(traits[which(traits$Species %in% s_i), "DP"], na.rm=T) / (sum(traits[which(traits$Species %in% s_i), "B"], na.rm=T)+1)
  indicators2[i,"Exo"] <- sum(traits[which(traits$Species %in% s_i), "NI"], na.rm=T)
  indicators2[i,"Commercial"] <- sum(traits[which(traits$Species %in% s_i), "all_commercial_level"], na.rm=T)
  indicators2[i,"High_commerc"] <- sum(traits[which(traits$Species %in% s_i), "highly_commercial_only"], na.rm=T)
}

###########################################################
## 4 Cryptobenthic (definition Brandl et al. 2018)
# Brandl SJ, Goatley CHR, Bellwood DR, Tornabene L. 2018 The hidden half: ecology and evolution of cryptobenthic fishes on coral reefs. Biol. Rev. 93, . (doi:10.1111/brv.12423))
###########################################################
crypto_families = c("Tripterygiidae", "Grammatidae", "Creediidae", "Aploactinidae", "Gobiidae", 
                    "Chaenopsidae", "Gobiesocidae", "Labrisomidae", "Pseudochromidae", "Bythitidae", 
                    "Plesiopidae", "Blenniidae", "Apogonidae", "Callionymidae", "Opistognathidae", "Syngnathidae")

crypto <- taxo %>%
  filter(Family %in% crypto_families) %>%
  pull(scientific_name)


for (i in 1:nrow(indicators2)) { # for each sample
  # list species present in the sample
  s_i <- data2 %>%
    filter(data2[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  crypto_i <- s_i[which(s_i %in% crypto)]
  
  # Calculate indicator
  indicators2[i,"Crypto"] <- length(crypto_i)
}

###########################################################
## 9 Chondrichtyens
###########################################################
for (i in 1:nrow(indicators2)) { # for each sample
  # list species present in the sample
  s_i <- data2 %>%
    filter(data2[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  elasmo_i <- s_i[which(s_i %in% elasmo)]
  
  # Calculate indicator
  indicators2[i,"Chondri"] <-length(elasmo_i)
}

##########################################################
## 6 - Community Temperature Index - CTI
## 10 - Mean trophic level
## 14 - Vulnerability
##########################################################
for (i in 1:nrow(indicators2)) { # for each sample
  # list species present in the sample
  s_i <- data2 %>%
    filter(data2[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # calculate indicators2
  indicators2[i,"CTI"] <- mean(traits[which(traits$Species %in% s_i), "SSTmean"], na.rm = T) 
  indicators2[i,"Trophic"] <- mean(traits[which(traits$Species %in% s_i), "Trophic_level"], na.rm=T) 
  indicators2[i,"Vulner"] <- mean(traits[which(traits$Species %in% s_i), "Vulnerability"], na.rm=T)
}

##########################################################
## 8 - Red List IUCN
##########################################################
### Make a dummy variable for IUCN categories
traits <- dummy_cols(traits, select_columns = 'IUCN_Red_List_Category')

## Calculate indicator
for (i in 1:nrow(indicators2)) { # for each sample
  # list species present in the sample
  s_i <- data2 %>%
    filter(data2[,i] == 1) %>%
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Number of species per category
  VU <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_VU"], na.rm=T)
  EN <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_EN"], na.rm=T)
  CR <- sum(traits[which(traits$Species %in% s_i), "IUCN_Red_List_Category_CR"], na.rm=T)
  
  # Calculate indicator
  indicators2[i,"RedList"] <- VU + EN + CR
}

##########################################################
## 13 - Phylogenetic Diversity - PD
##########################################################
# Retrieve the phylogeny of only native reef species across all three oceans.
phy <- fishtree_phylogeny(species = species)

plot(phy, show.tip.label = FALSE)
tiplabels(tip = which(phy$tip.label %in% species),
          pch=19, cex=2)

rownames(data2) <- gsub(" ", "_", rownames(data2), fixed = TRUE)

## check that phylogeny and data2 have matching names
#nc <- geiger::name.check(phy, data2) # 35 species not in tree
#missing <- gsub("_", " ", nc$data2_not_tree, fixed = TRUE)
## List of missing species that are not Chondrichtyens (as they are not in the tree)
#missing_fish <- traits[which(traits$Species %in% missing & traits$Super_class != "chondrichtyen"),"Species"]

# Manually check synonyms and find the species
species[species == "Mullus barbatus"] <- "Mullus barbatus barbatus"
rownames(data2)[rownames(data2) == "Mullus_barbatus"] <- "Mullus_barbatus_barbatus"

species[species == "Diplodus sargus"] <- "Diplodus sargus sargus"
rownames(data2)[rownames(data2) == "Diplodus_sargus"] <- "Diplodus_sargus_sargus"

species[species == "Diplodus cervinus"] <- "Diplodus cervinus cervinus"
rownames(data2)[rownames(data2) == "Diplodus_cervinus"] <- "Diplodus_cervinus_cervinus"

species[species == "Chelon auratus"] <- "Liza aurata"
rownames(data2)[rownames(data2) == "Chelon_auratus"] <- "Liza_aurata"

species[species == "Chelon ramada"] <- "Liza ramada"
rownames(data2)[rownames(data2) == "Chelon_ramada"] <- "Liza_ramada"

# Retrieve the missing phylogeny 
phy <- fishtree_phylogeny(species = species, type="phylogram")
nc <- geiger::name.check(phy, data2) # 40 espèces manquantes

# Remove from the data2 the species that are not in the tree
data22 <- data2[which(rownames(data2) %in% nc$data2_not_tree == F),]

# Transpose the ADNe matrix 
data22 <- t(data22)

# prune the tree
prunedTree <- prune.sample(data22,phy)

# Calculate PD
pd.result <- pd(data22, prunedTree, include.root=T)

# Add PD to indicator data2frame
indicators2[,"PD"] <- pd.result$PD

write.table(indicators2, file="01_Analyses_teleo/00_data/indicators_ports_2022_per_port.csv", sep=",")


#### Summarize total
indicators2  %>% 
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
png("01_Analyses_teleo/03_Outputs/boxplots_indicateurs_habitat.png", width = 900, height = 900)
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

