library(tidyverse)


# load csv spygen -> really messy !
data_metazoa <- read.csv("02_Analyses_metazoa/00_data/Metazoaires_SC21250_résultats_campagnes 1 et 2.csv", sep=";", na.strings = "")

# change column names 
colnames(data_metazoa) <- c(as.character(data_metazoa[4,1:6]), as.character(data_metazoa[3,7:ncol(data_metazoa)]))
colnames(data_metazoa)<- gsub(" ", "", colnames(data_metazoa))

# remove first rows
data_metazoa <- data_metazoa[-c(1:4),]

# keep only scientific_name and samples
data_metazoa <- data_metazoa[, -c(1:4,6)]

# replace empty cells with 0
for (i in 2:ncol(data_metazoa)) {
  data_metazoa[,i] <- gsub(" ", "", data_metazoa[,i])
}
data_metazoa[,2:ncol(data_metazoa)] <- as.numeric(unlist(data_metazoa[,2:ncol(data_metazoa)]))

data_metazoa[is.na(data_metazoa)] <- 0

# Remove rows with 0 obs
data_metazoa <- data_metazoa[rowSums(data_metazoa[,-1])!=0,]

# clean species_names
data_metazoa$scientific_name <- gsub(" ", "_", data_metazoa$scientific_name)

#save data
write.csv(data_metazoa, file="02_Analyses_metazoa/00_data/metazoa_reads.csv", row.names = F)

# Transform to presence_absence
data_metazoa_PA <- data_metazoa

data_metazoa_PA <- data_metazoa_PA %>% 
  mutate_if(is.numeric, ~1 * (. > 0))

# save data
write.csv(data_metazoa_PA, file="02_Analyses_metazoa/00_data/metazoa_presence.csv", row.names = F)
