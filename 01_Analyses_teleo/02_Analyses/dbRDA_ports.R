library(dplyr)
library(forcats)
library(stringr)
library(rsq)
library(margins)
library(betapart)
library(reshape)
library(tidyverse)
library(tidyselect)
library(vegan)
library(ggplot2)
library(patchwork)
library(ggalt)
library(ggrepel)
library(grid)
library(ggpubr)
library(RColorBrewer)
library(rdacca.hp)

###########################################################################################
## dbRDA ports only 
###########################################################################################
biodiv_port=read.csv("01_Analyses_teleo/00_data/matrice_teleo_port.csv", row.names=1) %>%
  t(.) %>%
  as.data.frame(.) %>%
  rownames_to_column(var='code_spygen')

# select species to keep
species_to_keep <- readRDS("01_Analyses_teleo/00_data/species_to_keep_ports.RDS")

biodiv_port <- biodiv_port[, c("code_spygen", species_to_keep)]

# metadata
meta_port=read.csv("00_Metadata/metadata_port.csv", sep=";") 

# combiner toutes les donnees
data_port=left_join(biodiv_port,meta_port)
rownames(data_port)=data_port[,"code_spygen"]


# remove biohut
data_port = data_port %>%
  filter(type == "Port")

#### faire la dbRDA
data_dbrda_port <- data_port[,c(2:123)]
meta_dbrda_port <- data_port[,c(124:ncol(data_port))]


# Hierarchical and Variation Partitioning for Canonical Analysis 
# Include all variable - Table Sxx


rdacca.hp(vegdist(data_dbrda_port,method="jaccard"),
          meta_dbrda_port[,c("Campaign","Certification","Area_ha","Depth_m","Habitat", "Longitude")],method="dbRDA",add=F,type="R2")


rda.port<-rdacca.hp(vegdist(data_dbrda_port,method="jaccard"),
                    meta_dbrda_port[,c("Campaign","Certification","Area_ha","Depth_m","Habitat", "Longitude")],method="dbRDA",add=F,type="R2")

rda.port
plot(rda.port)
ggsave(file="01_Analyses_teleo/03_Outputs/rdacca_port_FigureSXX.png")

rda_port<- as.data.frame(rda.port$Hier.part)
write.csv(rda_port, file="01_Analyses_teleo/03_Outputs/Table_SXX_rdacca_port.csv", row.names = F)


permu.hp(vegdist(data_dbrda_port,method="jaccard"),
         meta_dbrda_port[,c("Campaign","Certification","Area_ha","Depth_m","Habitat", "Longitude")],method="dbRDA",add=F,type="R2",permutations=100)

permu.hp(vegdist(data_dbrda_port,method="jaccard"),
         meta_dbrda_port[,c("Campaign","Certification","Area_ha","Depth_m","Habitat", "Longitude")],method="dbRDA",add=F,type="R2",permutations=1000)



#### dbRDA Condition(longitude)

RDA_port<- capscale(data_dbrda_port ~ Certification + Area_ha + Depth_m + Habitat + Condition(Longitude), 
                    meta_dbrda_port, distance ="jaccard")

# get scores
site_scores <- scores(RDA_port)$sites     ## separating out the site scores
species_scores <- scores(RDA_port)$species %>% data.frame()   ## separating out the species
rsqr <- round(RsquareAdj(RDA_port)$adj.r.squared*100, 2) # r squared

# extract the percentage variability explained by axes
sumdbrda <- summary(RDA_port)

### Check the collinearity of the model
vif.cca(RDA_port)

### Test the significance
anova(RDA_port, perm=999)
anova(RDA_port, perm=999, by="margin") # significant predictors : Area, Depth, Habitat
anova(RDA_port, perm=999, by="axis")

### Getting the scores for plotting.
scrs <- scores(RDA_port, display = c("species", "sites", "bp"), choices = 1:4, scaling = 2)
scrs

### Collect information on the pcao axes
species_centroids <- data.frame(scrs$species)
species_centroids
species_centroids$PC_names <- rownames(species_centroids)

### Add the PC components names.
species_centroids$PC_names <- rownames(species_centroids) 

### Add information of arrows
continuous_arrows <- data.frame(scrs$biplot)
continuous_arrows
rownames(continuous_arrows) <- c("Certification", "Area", "Depth", "Habitat")

### Add names of variables
continuous_arrows$number <- rownames(continuous_arrows)
continuous_arrows
baseplot<-plot(RDA_port, scaling = 2)
mult <- attributes(baseplot$biplot)$arrow.mul

### Add port names
sitenames<- data_port %>%
  arrange(match(code_spygen, rownames(scrs$sites))) %>%
  pull(site)
sites_centroids <- data.frame(scrs$sites)
sites_centroids$SITE <- sitenames
head(sites_centroids)

# order ports from east to west
sites_centroids$SITE <- factor(sites_centroids$SITE,                                   
                               levels = c("Cannes","Porquerolles", "La Ciotat",
                                          "Saintes Maries de la Mer", "Marseillan",
                                          "Agde", "Port Vendres"))

### Make a ggplot
RDA_plot <- ggplot(data = sites_centroids, aes(x = CAP1, y = CAP2))+
  geom_hline(yintercept = 0, lty = "dotted") + geom_vline(xintercept = 0, lty = "dotted") +
  geom_point(data = sites_centroids, pch = 21, size = 4, aes(fill = SITE))+
  scale_fill_manual(values = brewer.pal(7,"RdYlBu") ,name="Port - East to West")+
  xlim(-2,2.3) +
  geom_segment(data = continuous_arrows,
               aes(x = 0, xend = mult * CAP1, y = 0, yend = mult * CAP2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = continuous_arrows,
            aes(x= (mult + mult/10) * CAP1, y = (mult + mult/10) * CAP2,
                label = number), size = 4, hjust = 0.5, vjust = -0.2)+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = paste("RDA1 (", round(sumdbrda$cont$importance[2,1]*100, 2), "%)", sep = ""), y = paste("RDA2 (", round(sumdbrda$cont$importance[2,2]*100, 2), "%)", sep = ""))+
  theme(axis.text = element_text(colour = "black", size = 12)) +
  theme(axis.title = element_text(size = 12, colour = "black", family = "Helvetica", face = "bold")) +
  theme_classic()
RDA_plot

# export figure
ggsave(plot = RDA_plot, filename = "01_Analyses_teleo/03_Outputs/dbRDA_ports_partiel_longitude.jpeg", 
       dpi = 600)
