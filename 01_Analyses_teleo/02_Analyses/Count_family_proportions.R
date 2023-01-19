library(tidyverse)
library(ggplot2)
library(rfishbase)
library(reshape)
'%ni%' <- Negate("%in%")

# load fish taxonomy
x <- load_taxa()
taxo <- collect(x)


# load data port
species_port <- read.csv("01_Analyses_teleo/00_data/matrice_teleo_port.csv", row.names=1) %>%
  row.names(.) %>%
  as.data.frame(.)

colnames(species_port) <- "species"
species_port$port <- 1


# load outside data

species_outside <- read.csv("01_Analyses_teleo/00_data/biodiv_milieu_naturel.csv", row.names = 1) %>%
  row.names(.) %>%
  as.data.frame(.)

colnames(species_outside) <- "species"
species_outside$outside <- 1
species_outside$species <- gsub(" ", "_", species_outside$species)


# join both dataset

species_all <- full_join(species_port, species_outside)
species_all[is.na(species_all)] <- 0

for (i in 1:nrow(species_all)) {
  if (species_all[i,"port"]==1 & species_all[i,"outside"]==1) {
    species_all[i,"both"] <- 1
    species_all[i,"port"] <- 0
    species_all[i,"outside"] <- 0
    }
}
species_all[is.na(species_all)] <- 0


# filter genus assignments and add family
genus <- species_all %>% 
  filter(grepl("_sp.$", species))

species_all <- species_all %>%
  filter(species %ni% genus$species)

genus$genus <- sub("_sp.", "", genus$species)

genus <- left_join(genus, taxo[,c("Genus", "Family")], by=c("genus"="Genus"))%>%
  distinct(species, .keep_all=T)

genus <- genus[,c(6,2,3,4)] 


# filter species assignments and add family
species <- species_all %>% 
  filter(grepl("_", species))

species$species <- gsub("_", " ", species$species)

species <- left_join(species, taxo[,c("Species", "Family")], by=c("species"="Species"))%>%
  distinct(species, .keep_all=T)

species_NA <- species %>% 
  filter(is.na(Family))

species <- species %>% filter(species %ni% species_NA$species)

species_NA <- species_NA %>%
  mutate(Family = case_when(
    species=="Parablennius incognitus P sanguinolentus" ~ "Blenniidae",
    species=="Labrus merula Labrus viridis" ~ "Labridae",
    species=="C. lucerna L. dieuzeidei" ~ "Triglidae",
    species=="C. obscurus T. lastoviza" ~ "Triglidae",
    species=="E. gurnardus T. lyra" ~ "Triglidae",
    species=="Spicara flexuosa S smaris" ~ "Sparidae",
    species=="B. boops O. melanura" ~ "Labridae",
    species=="D. dentex P. auriga P. pagrus" ~ "Labridae",
    species=="Hexanchus griseus" ~ "Hexanchidae",
    species=="Raja asterias Raja clavata" ~ "Rajidae",
    species=="C. heterurus H. speculiger" ~ "Exocoetidae",
    species=="Tripterygion delaisi xanthosoma" ~ "Tripterygiidae",
    species=="Notoscopelus elongatus kroyeri" ~ "Myctophidae",
    species=="Spicara flexuosa S maena" ~ "Sparidae",
    species=="Isurus oxyrinchus" ~ "Lamnidae"
  ))

species <- rbind(species, species_NA)
species <- species[,c(5,2,3,4)] 

# keep family assignments

family <- species_all %>% 
  filter(!grepl("_", species))

colnames(family) <- c("Family", "port", "outside", "both")

# Bind all datasets and count number of family occurrences

family_all <- rbind(family, genus, species)

fam <- unique(family_all$Family)
fam_proportion <- data.frame(Family=character(102), port=numeric(102), outside=numeric(102), both=numeric(102))

for (i in 1:length(fam)) {
  df <- family_all %>%
    filter(Family==fam[i])
  fam_proportion[i,"Family"] <- fam[i]
  fam_proportion[i,"port"] <- sum(df$port)
  fam_proportion[i,"outside"] <- sum(df$outside)
  fam_proportion[i,"both"] <- sum(df$both)
}

write.csv(fam_proportion, file="01_Analyses_teleo/03_Outputs/count_families.csv", row.names = F)


#number of families in ports
names(fam_proportion)
n_port = fam_proportion[fam_proportion$port>0 | fam_proportion$both>0, ]
dim(n_port) #54 families


# Plot graph count families

fam_proportion2 <- melt(fam_proportion)
colnames(fam_proportion2) <- c("Family", "Zone", "Count")


plot <- ggplot(fam_proportion2, aes(x=reorder(Family, Count), y = Count, fill = Zone)) + 
  geom_bar(stat="identity", show.legend = TRUE) + 
  theme_bw() +
  scale_fill_manual(values =c("#FFC47E","#7FC5C5", "grey"))+ 
  labs(x="Family", y="Number of species")+
  coord_flip()


ggsave(plot, file="01_Analyses_teleo/03_Outputs/Families_count.png", width=6, height=16)
