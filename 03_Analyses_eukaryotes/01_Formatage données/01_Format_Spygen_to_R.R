library(tidyverse)

###################################################################################################################
# Eukaryote campaign 1
###################################################################################################################

# load csv spygen -> really messy !
data_euk1 <- read.csv("03_Analyses_eukaryotes/00_data/Euka_SC21250_CEFE_campagne 1.csv", sep=";", na.strings = "")

# change column names 
colnames(data_euk1) <- c(colnames(data_euk1[,1:4]), as.character(data_euk1[1,5:ncol(data_euk1)]))

# remove first row
data_euk1 <- data_euk1[-1,]

# keep only class_order
data_euk1$order <- gsub("-",NA,data_euk1$order)
data_euk1$class <- gsub("-",NA,data_euk1$class)
data_euk1$class_order <- paste(data_euk1$class, data_euk1$order, sep="_")
data_euk1 <- data_euk1 %>%
  filter(class_order!="NA_NA")

data_euk1 <- data_euk1 %>%
  select(class_order, c(5:19))

# replace empty cells with 0
data_euk1[is.na(data_euk1)] <- 0

# Remove rows with 0 obs
data_euk1[,2:ncol(data_euk1)] <- as.numeric(unlist(data_euk1[,2:ncol(data_euk1)]))
data_euk1 <- data_euk1[rowSums(data_euk1[,-1])!=0,]


###################################################################################################################
# Eukaryote campaign 2
###################################################################################################################


# load csv spygen -> really messy !
data_euk2 <- read.csv("03_Analyses_eukaryotes/00_data/Euka_SC21250_CEFE_campagne 2.csv", sep=";", na.strings = "")

# change column names 
colnames(data_euk2) <- as.character(data_euk2[1,])

# remove first row
data_euk2 <- data_euk2[-1,]

# keep only class_order
data_euk2$order <- gsub("-",NA,data_euk2$order)
data_euk2$class <- gsub("-",NA,data_euk2$class)
data_euk2$class_order <- paste(data_euk2$class, data_euk2$order, sep="_")
data_euk2 <- data_euk2 %>%
  filter(class_order!="NA_NA")

data_euk2 <- data_euk2 %>%
  select(class_order, c(5:18))

# replace empty cells with 0
data_euk2[is.na(data_euk2)] <- 0

# Remove rows with 0 obs
data_euk2[,2:ncol(data_euk2)] <- as.numeric(unlist(data_euk2[,2:ncol(data_euk2)]))
data_euk2[is.na(data_euk2)] <- 0
data_euk2 <- data_euk2[rowSums(data_euk2[,-1])!=0,]

###########################################################################################################
# Assemble two data frames

data_euk <- full_join(data_euk1, data_euk2)

#save data
write.csv(data_euk, file="03_Analyses_eukaryotes/00_data/eukaryote_reads.csv", row.names = F)

# Transform to presence_absence
data_euk_PA <- data_euk

data_euk_PA <- data_euk_PA %>% 
  mutate_if(is.numeric, ~1 * (. > 0))

# save data
write.csv(data_euk_PA, file="03_Analyses_eukaryotes/00_data/eukaryote_presence.csv", row.names = F)
