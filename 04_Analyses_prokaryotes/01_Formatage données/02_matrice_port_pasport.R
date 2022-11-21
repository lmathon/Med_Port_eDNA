### Load the libraries
library(dplyr)
library(rsq)
library(purrr)
library(ggplot2)
library(vegan)
library(betapart)
library(textshape)


### Load and clean ADNe eukaryote hors ports
 
adne <- read.csv("04_Analyses_prokaryotes/00_data/Prokaryotes_milieu_naturel.csv", header=T, sep=";", na.strings = "")
colnames(adne) <- c(colnames(adne[,1:6]), as.character(adne[1,7:ncol(adne)]))
adne <- adne[-c(1, nrow(adne)),]

# add column class_order
adne <- adne %>%
  select(c(6:ncol(adne)))
adne$Taxon <- gsub(" ", "_", adne$Taxon)
colnames(adne)[1] <- "taxon"
adne <- adne %>%
  filter(!is.na(taxon))

adne[,2:ncol(adne)] <- as.numeric(unlist(adne[,2:ncol(adne)]))
adne[is.na(adne)] <- 0

# keep samples for which we have metadata
meta_nat <- read.csv("00_Metadata/metadata_milieu_naturel.csv")
keep <- c("taxon", as.character(meta_nat$SPYGEN_code))

adne <- adne[, which(names(adne) %in% keep)]


###ADNe port
adne_port <- read.csv("04_Analyses_prokaryotes/00_data/prokaryote_presence.csv", header = T) #--> prends les noms de ligne pour en faire une colonne 

#jointure by Species 
adne_tot <- adne %>%
  full_join(adne_port) %>% #NA quand sp. pas dans les matrices 
  replace(is.na(.), 0)

adne_tot <- adne_tot %>% 
  mutate_if(is.numeric, ~1 * (. > 0))


# Species for rownames
rownames(adne_port) <- adne_port$taxon
adne_port <- adne_port[,-1]

write.csv(adne_port, "04_Analyses_prokaryotes/00_data/matrice_prokaryote_port.csv")



rownames(adne_tot) <- adne_tot$class_order
adne_tot <- adne_tot[,-1]

write.csv(adne_tot, "04_Analyses_prokaryotes/00_data/matrice_prokaryote_totale.csv") 




