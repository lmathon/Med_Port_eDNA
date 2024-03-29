### Load the libraries
library(dplyr)
library(rsq)
library(purrr)
library(ggplot2)
library(vegan)
library(betapart)
library(textshape)

###ADNe
adne_port <- read.csv("01_Analyses_teleo/00_data/teleo_presence.csv", header = T) #--> prends les noms de ligne pour en faire une colonne 
colnames(adne_port)[1] <- "Species"

 

adne <- read.csv("01_Analyses_teleo/00_data/biodiv_milieu_naturel.csv", header=T, row.names=1) %>%
  tibble::rownames_to_column(var="Species") #--> prends les noms de ligne pour en faire une colonne 
adne$Species <- gsub(" ", "_", adne$Species)



# Species for rownames
species_to_keep <- readRDS("01_Analyses_teleo/00_data/species_to_keep_ports.RDS")

adne_port <- adne_port %>%
  filter(Species %in% species_to_keep)

rownames(adne_port) <- adne_port$Species
adne_port <- adne_port[,-1]

write.csv(adne_port, "01_Analyses_teleo/00_data/matrice_teleo_port.csv")

#jointure by Species 
adne_port$Species <- rownames(adne_port)

adne_tot <- adne %>%
  full_join(adne_port, by = "Species") %>% #NA quand sp. pas dans les matrices 
  replace(is.na(.), 0)


rownames(adne_tot) <- adne_tot$Species
adne_tot <- adne_tot[,-1]

write.csv(adne_tot, "01_Analyses_teleo/00_data/matrice_teleo_totale.csv") 




