library(tidyverse)
library(stringr)


### METADATA Port

# load October 2021

meta_port <- read.csv("00_Metadata/metadata_eDNA_port_october_2021.csv", header=T)

meta_port <- meta_port[,1:9]
meta_port <- meta_port[,-c(2,3,7)]

meta_port$port <- "true"

# load June 2022

meta_port2 <- read.csv("00_Metadata/metadata_eDNA_port_june_2022.csv", sep=";", header=T)

meta_port2 <- meta_port2[1:14,1:9]
meta_port2 <- meta_port2[,-c(2,3,7)]

meta_port2$port <- "true"

# Assemble 2 seasons
meta_port <- rbind(meta_port, meta_port2)

# add certification "port propre" for Cannes, La Ciotat, Sainte-Marie-de-la-mer and Porquerolles
meta_port$port_propre <- NA

meta_port <- meta_port %>%
  mutate(port_propre = case_when(
    site=="Cannes" ~ "Port propre",
    site=="La Ciotat" ~ "Port propre",
    site=="Porquerolles" ~ "Port propre",
    site=="Saintes Maries de la Mer" ~ "Port propre",
    site=="Marseillan " ~ "Non certifie",
    site=="Agde" ~ "Non certifie",
    site=="Port Vendres" ~ "Non certifie"
    ))

meta_port$site <- str_trim(meta_port$site, "right")
write.csv(meta_port, "00_Metadata/metadata_port.csv", row.names = F) 



### METADATA Outside

meta2 <- read.csv("00_Metadata/metadata_milieu_naturel.csv", header=T)


meta2 <- meta2[,c('SPYGEN_code','Date', 'Site', 'Site', 'projet', 'protection')]
colnames(meta2) <- c('code_spygen','date','station','site','project','habitat')
meta2$port <-"false"

meta2 <- meta2 %>%
  mutate(site = case_when(
    station=="3" ~ "Banyuls_outside",
    station=="5" ~ "Banyuls_outside",
    station=="1" ~ "Banyuls_outside",
    station=="2" ~ "Banyuls_outside",
    station=="4" ~ "Banyuls_outside",
    station=="carry-le-Rouet" ~ "carry-le-Rouet_outside",
    station=="porquerolles" ~ "porquerolles_outside",
    station=="cap_roux" ~ "cap-roux_outside",
    station=="cap_de_nice" ~ "cap-de-nice_outside",
    station=="baie _de_villefranche" ~ "baie-de-villefranche_outside",
    station=="Calanques" ~ "Calanques_outside",
    station=="beach_rock_cassis" ~ "Cassis_outside",
    station=="cap_lardier" ~ "cap-lardier_outside",
    station=="cap_negre" ~ "cap-negre_outside",
    station=="Banyuls" ~ "Banyuls_outside",
    station=="banyuls" ~ "Banyuls_outside",
    station=="embiez" ~ "Embiez_outside",
    station=="carry" ~ "carry-le-Rouet_outside",
    station=="calanques" ~ "Calanques_outside",
    station=="Calvi" ~ "Calvi_outside"
      ))

meta_port <- meta_port %>%
  mutate(site = case_when(
    site=="Cannes" ~ "Cannes_port",
    site=="La Ciotat" ~ "La Ciotat_port",
    site=="Marseillan" ~ "Marseillan_port",
    site=="Saintes Maries de la Mer" ~ "Saintes Maries de la Mer_port",
    site=="Agde" ~ "Agde_port",
    site=="Port Vendres" ~ "Port Vendres_port",
    site=="Porquerolles" ~ "Porquerolles_port"
  ))



meta_tot <- rbind(meta2,meta_port[,-8])

write.csv(meta_tot, "00_Metadata/metadata_tot.csv", row.names = F) 


