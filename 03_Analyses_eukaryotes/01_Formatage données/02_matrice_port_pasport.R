### Load the libraries
library(dplyr)
library(rsq)
library(purrr)
library(ggplot2)
library(vegan)
library(betapart)
library(textshape)


### Load and clean ADNe eukaryote hors ports
 
adne <- read.csv("03_Analyses_eukaryotes/00_data/Eukaryotes_milieu_naturel.csv", header=T, sep=";", na.strings = "")
colnames(adne) <- c(colnames(adne[,1:4]), as.character(adne[1,5:ncol(adne)]))
adne <- adne[-c(1, nrow(adne)),]

# add column class_order
adne$class_order <- paste(adne$class, adne$order, sep="_")
adne <- adne %>%
  filter(class_order!="NA_NA")

adne <- adne %>%
  select(class_order, c(5:16))

adne[,2:ncol(adne)] <- as.numeric(unlist(adne[,2:ncol(adne)]))
adne[is.na(adne)] <- 0

# keep samples for which we have metadata
meta_nat <- read.csv("00_Metadata/metadata_milieu_naturel.csv")
keep <- c("class_order", as.character(meta_nat$SPYGEN_code))

adne <- adne[, which(names(adne) %in% keep)]

# deal with duplicated identifications
dup1 <- adne %>%
  filter(class_order=="Dinophyceae_Peridiniales")

dup1$class_order <- c("Dinophyceae_Peridiniales1", "Dinophyceae_Peridiniales2")

adne <- adne %>%
  filter(class_order!="Dinophyceae_Peridiniales")
adne <- rbind(adne, dup1)


dup2 <- adne %>%
  filter(class_order=="NA_Thaumatomonadida")

dup2$class_order <- c("NA_Thaumatomonadida1", "NA_Thaumatomonadida2")

adne <- adne %>%
  filter(class_order!="NA_Thaumatomonadida")
adne <- rbind(adne, dup2)

###ADNe port
adne_port <- read.csv("03_Analyses_eukaryotes/00_data/eukaryote_presence.csv", header = T) #--> prends les noms de ligne pour en faire une colonne 

#jointure by Species 
adne_tot <- adne %>%
  full_join(adne_port) %>% #NA quand sp. pas dans les matrices 
  replace(is.na(.), 0)

adne_tot <- adne_tot %>% 
  mutate_if(is.numeric, ~1 * (. > 0))


# Species for rownames
rownames(adne_port) <- adne_port$class_order
adne_port <- adne_port[,-1]

write.csv(adne_port, "03_Analyses_eukaryotes/00_data/matrice_eukaryote_port.csv")



rownames(adne_tot) <- adne_tot$class_order
adne_tot <- adne_tot[,-1]

write.csv(adne_tot, "03_Analyses_eukaryotes/00_data/matrice_eukaryote_totale.csv") 




