library(tidyverse)

###################################################################################################################
# prokaaryote campaign 1
###################################################################################################################

# load csv spygen -> really messy !
data_proka1 <- read.csv("04_Analyses_prokaryotes/00_data/Bact_SC21250_CEFE_campagne 1.csv", sep=";", na.strings = "")

# change column names 
colnames(data_proka1) <- c(colnames(data_proka1[,1:5]), as.character(data_proka1[1,6:ncol(data_proka1)]))

# remove first row
data_proka1 <- data_proka1[-1,]

# keep only taxon
data_proka1 <- data_proka1 %>%
  select(c(5:ncol(data_proka1)))
data_proka1$taxon <- gsub(" ", "_", data_proka1$taxon)
data_proka1 <- data_proka1 %>%
  filter(!is.na(taxon))

# replace empty cells with 0
data_proka1[,2:ncol(data_proka1)] <- as.numeric(unlist(data_proka1[,2:ncol(data_proka1)]))
data_proka1[is.na(data_proka1)] <- 0

# Remove rows with 0 obs
data_proka1 <- data_proka1[rowSums(data_proka1[,-1])!=0,]


###################################################################################################################
# prokaaryote campaign 2
###################################################################################################################


# load csv spygen -> really messy !
data_proka2 <- read.csv("04_Analyses_prokaryotes/00_data/Bact_SC21250_CEFE_campagne 2.csv", sep=";", na.strings = "")

# change column names 
colnames(data_proka2) <- as.character(data_proka2[1,])
colnames(data_proka2) <- gsub(" ", "", colnames(data_proka2))

# remove first row
data_proka2 <- data_proka2[-1,]

# keep only scientific_name
data_proka2 <- data_proka2 %>%
  select(c(5:ncol(data_proka2)))
data_proka2$scientific_name <- gsub(" ", "_", data_proka2$scientific_name)
colnames(data_proka2)[1] <- "taxon"
data_proka2 <- data_proka2 %>%
  filter(!is.na(taxon))

# replace empty cells with 0
for (i in 2:ncol(data_proka2)) {
  data_proka2[,i] <- gsub(" ", "", data_proka2[,i])
}

data_proka2[,2:ncol(data_proka2)] <- as.numeric(unlist(data_proka2[,2:ncol(data_proka2)]))
data_proka2[is.na(data_proka2)] <- 0

# Remove rows with 0 obs
data_proka2 <- data_proka2[rowSums(data_proka2[,-1])!=0,]

###########################################################################################################
# Assemble two data frames

data_proka <- full_join(data_proka1, data_proka2)
data_proka[is.na(data_proka)] <- 0


#save data
write.csv(data_proka, file="04_Analyses_prokaryotes/00_data/prokaryote_reads.csv", row.names = F)

# Transform to presence_absence
data_proka_PA <- data_proka

data_proka_PA <- data_proka_PA %>% 
  mutate_if(is.numeric, ~1 * (. > 0))

# save data
write.csv(data_proka_PA, file="04_Analyses_prokaryotes/00_data/prokaryote_presence.csv", row.names = F)
