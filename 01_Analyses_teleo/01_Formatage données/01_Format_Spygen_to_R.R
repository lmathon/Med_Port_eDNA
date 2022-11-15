library(tidyverse)


# load csv spygen -> really messy !
data_teleo <- read.csv("01_Analyses_teleo/00_data/Teleo_SC21250_résultats_campagnes 1 et 2.csv", sep=";")

# Cas particulier : code spygen manquant dans la colonne "nb seq" du premier echantillon
data_teleo[3,7] <- data_teleo[3,6]

# Save the taxo in a new DF
taxo <- data_teleo[5:nrow(data_teleo), 1:4]
colnames(taxo) <- c("Class", "Order", "Family", "scientific_name")
#save taxo
write.csv(taxo, file="01_Analyses_teleo/00_data/teleo_taxo.csv", row.names = F)

# change column names 
colnames(data_teleo) <- c(as.character(data_teleo[4,1:5]), as.character(data_teleo[3,6:ncol(data_teleo)]))
colnames(data_teleo)<- gsub(" ", "", colnames(data_teleo))

# remove columns "Nbr rep"
cols_to_remove <- numeric()

for (i in 1:ncol(data_teleo)) {
  if (data_teleo[4,i]==" Nbr rep "){
    cols_to_remove <- c(cols_to_remove, i)
  }
}

data_teleo <- data_teleo %>%
  select(-cols_to_remove)


# remove first rows
data_teleo <- data_teleo[-c(1:4),]

# keep only scientific_name and samples
data_teleo <- data_teleo[, -c(1:3,5)]
data_teleo[,2:ncol(data_teleo)] <- as.numeric(unlist(data_teleo[,2:ncol(data_teleo)]))


# replace empty cells with 0
data_teleo[is.na(data_teleo)] <- 0

# Remove rows with 0 obs
data_teleo <- data_teleo[rowSums(data_teleo[,-1])!=0,]

# clean species_names
data_teleo$scientific_name <- gsub(" ", "_", data_teleo$scientific_name)

# Resolve duplicated species names (Engraulis_encrasicolus)
dup <- data_teleo %>%
  filter(scientific_name=="Engraulis_encrasicolus")
new <- data.frame(t(c("Engraulis_encrasicolus", colSums(dup[,-1]))))
colnames(new)[1] <- "scientific_name"


data_teleo <- data_teleo %>%
  filter(scientific_name!="Engraulis_encrasicolus") %>%
  rbind(new)

data_teleo[,2:ncol(data_teleo)] <- as.numeric(unlist(data_teleo[,2:ncol(data_teleo)]))

#save data
write.csv(data_teleo, file="01_Analyses_teleo/00_data/teleo_reads.csv", row.names = F)

# Transform to presence_absence
data_teleo_PA <- data_teleo

data_teleo_PA <- data_teleo_PA %>% 
  mutate_if(is.numeric, ~1 * (. > 0))

# save data
write.csv(data_teleo_PA, file="01_Analyses_teleo/00_data/teleo_presence.csv", row.names = F)
