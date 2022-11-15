### Load the libraries
library(dplyr)
library(rsq)
library(purrr)
library(ggplot2)
library(vegan)
library(betapart)


### Load Metazoaires port
adne_port <- read.csv("02_Analyses_metazoa/00_data/metazoa_presence.csv", header = T, row.names=1) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var="Species") #--> prends les noms de ligne pour en faire une colonne 
adne_port$Species <- sub(pattern = " ",  replacement = "_",  adne_port$Species)

write.csv(adne_port, "02_Analyses_metazoa/00_data/matrice_metazoa_port.csv", row.names = F) 



# Load Metazoaires outside 
adne <- read.csv("02_Analyses_metazoa/00_data/biodiv_milieu_naturel.csv", header=T, row.names=1) %>%
  t(.) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var="Species") #--> prends les noms de ligne pour en faire une colonne 

#jointure by Species 
adne_tot <- adne %>%
  full_join(adne_port, by = "Species") %>% #NA quand sp. pas dans les matrices 
  replace(is.na(.), 0) %>%
  column_to_rownames(var = "Species")

write.csv(adne_tot, "02_Analyses_metazoa/00_data/matrice_metazoa_totale.csv") 




