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
data <- read.csv("01_Analyses_teleo/00_data/teleo_presence.csv", header=T, row.names=1) %>%
  filter(rowSums(.) > 0) %>%
  t(.)  %>%
  as.data.frame(.)

meta <- read.csv("00_Metadata/metadata_port.csv", header=T)

traits <- read.csv("01_Analyses_teleo/00_data/Functional_data_corrected_20220124.csv", header=T)


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
loadNamespace("rfishbase") #installer manuellement le package rfishbase (via tools)


# Quelles espèces ne sont pas dans les traits ?
check <- species[which(species %in% traits$Species ==F)]
print(check)
#on en a trouvé 12 qui n'étaient pas dans la matrice des traits 
# Remove the species from the data for now 
data2 <- data %>%
  filter(rownames(.) %in% check ==F)
write.csv(data2, file="01_Analyses_teleo/00_data/teleo_presence_in_traits.csv")

#######################################################################################################
## Calculate the indicators for each sample
#######################################################################################################
## Create the result matrix
indicators <- matrix(NA,ncol(data), 14, #choix du nombre de colonnes de data pour le nombre de lignes de indicators, choix de 14 colonnes
                     dimnames=list(colnames(data), #noms des lignes = noms des colonnes de data 
                                   c("R",  "FD", "LFI", "Crypto", "DP_B_ratio",  "CTI", "Exo","RedList", "Chondri", "Trophic", "Commercial", "High_commerc", "PD", "Vulner"))) #choix manuel du nom des lignes

##########################################################
## 1 - Species Richness R
##########################################################
indicators[,1] <- apply(data, 2, sum) #dans la colonne 1 de indicator on met la somme du nb d'espèces, 2 = dimension de la matrice 

##########################################################
## 2 - Functional diversity FD
##########################################################
for (i in 1:nrow(indicators)) { # for each sample = chaque ligne de indicators --> boucle n°1
  # list species present in the sample
  s_i <- data %>%
    dplyr::filter(data[,i] == 1) %>% #pour toutes les sp présentes = valeur de 1
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp)
  
  # Get the functional groups of these species
  fd_i <- as.factor(traits[which(traits$Species %in% s_i), "m.cluster_core"]) #mettre sous forme de facteur les traits pour lesquels l'argument de la colonne species appartient à s_i, soit si l'espèces est bien présente (valeur 1)
  #cluster_core correspond aux groupes fonctionnels 
  
  # Number of unique functional groups
  indicators[i,2] <- length(unique(fd_i)) # diversité fonctionnelle = longueur du facteur fd_i qui contient les groupes fonctionnels sans compter les doublons
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
  indicators[i,"DP_B_ratio"] <- sum(traits[which(traits$Species %in% s_i), "DP"]) / (sum(traits[which(traits$Species %in% s_i), "B"])+1)
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
  mutate(crypto_Brandl = if_else(Family %in% crypto_families, 1,0)) #si family in crypto_families, vaut 1, sinon 0 ; transforme un dataframe existant 

for (i in 1:nrow(indicators)) { # for each sample = chaque ligne 
  # list species present in the sample
  s_i <- data %>%
    filter(data[,i] == 1) %>% #si sp bien présente
    tibble::rownames_to_column(var="Sp") %>%
    pull(Sp) #extraire une unique colonne
  
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
traits <- dummy_cols(traits, select_columns = 'IUCN_Red_List_Category') #dummy = binary, creation de colonnes binaires pour chaque valeur possible = catégorie IUCN 

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
phy <- fishtree_phylogeny(species = species) #fishtree phylogeny c'est une database publique a priori 
#Requested 112 but only found 82 species.

plot(phy, show.tip.label = FALSE)
tiplabels(tip = which(phy$tip.label %in% species),
          pch=19, cex=2) #ajoute des étiquettes aux noeuds; pch = type de symbole, cex = size

rownames(data) <- gsub(" ", "_", rownames(data), fixed = TRUE) #remplacement des matchs "à remplacer par", "ce symbole"

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
nc <- geiger::name.check(phy, data) # 30 espèces manquantes

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

write.csv(indicators, file="01_Analyses_teleo/00_data/indicators_ports.csv")

#### Summarize indicators per habitat (biohut ou non) per port
indicators_port_habitat = indicators  %>% 
  as.data.frame(.) %>%
  rownames_to_column(var="code_spygen") %>%
  inner_join(meta, by="code_spygen") %>%
  #selectionner les sites où les biohuts on été échantillonné
  filter(site %in% c("Cannes", "Marseillan ", "Saintes Maries de la Mer")) %>%
  select(R, FD, Crypto, Chondri, Commercial, Exo, RedList, habitat, site) %>%
  group_by(site, habitat) %>%
  
  summarise(across(everything(), list(mean))) #résumé parmi plusieurs colonnes





# Charger datas ports
meta <- read.csv("00_Metadata/metadata_port.csv", header=T)
ind_ports <- read.csv("01_Analyses_teleo/00_data/indicators_ports.csv", header=T, row.names=1) %>%
  mutate(Location = "port")

# Charger datas milieu naturel (= réserve et hors réserve)
meta_nat <- read.csv("00_Metadata/metadata_milieu_naturel.csv", header=T)
ind_nat <- read.csv("01_Analyses_teleo/00_data/indicators_milieu_naturel.csv", header=T, row.names=1) %>%
  rownames_to_column(var="spygen_code") %>%
  inner_join(meta_nat[,c("spygen_code", "protection")], by="spygen_code") %>%
  rename(Location = protection) %>%
  column_to_rownames(var="spygen_code") %>%
  select(colnames(ind_ports))

# Combiner les deux datasets
ind_all <- rbind(ind_nat,ind_ports) %>%
  mutate(DP_B_ratio = log10(DP_B_ratio))

# Boxplots
# Renommer les points milieu naturel = outside et reserve en "hors port"
ind_all$Location <- sub(pattern = "outside",  replacement = "hors port", ind_all$Location)
ind_all$Location <- sub(pattern = "reserve",  replacement = "hors port", ind_all$Location)

ind_all$Location <- factor(ind_all$Location,                                    # Factor levels in decreasing order
                            levels = c("hors port","port"))


p <-list()
## Draw the plot
l=1
for (i in c(1,2,8,11)) { # Choix des indicateurs à présenter
  # Gather data
  dat_i <- ind_all[,c(i,15)]
  colnames(dat_i)[1] <- "Y"
  
  p[[l]]<-ggplot(data=dat_i, aes(x=Location, y=Y, color=Location)) +
  geom_boxplot(notch=F) +
  scale_color_manual(values=c("darkblue", "cyan4",  "darkorange")) +
  ylab(colnames(ind_all)[i]) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.y=element_text(colour="black",size=16)) +
  theme(axis.text.x=element_text(colour="black",size=18)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_text(colour="black",size=18))+ 
  #stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size=4)
  stat_compare_means(aes(label = paste0("p = ", ..p.signif..)), size=6)

  l=l+1  
}

## Save plot
png("01_Analyses_teleo/03_Outputs/Boxplots_indicators_port_milieu_nat3.png", 
    width = 900, height = 1000) 
do.call(grid.arrange,c(p, list(ncol=2)))
dev.off()

