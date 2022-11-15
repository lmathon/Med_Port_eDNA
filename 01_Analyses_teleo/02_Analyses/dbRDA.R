library(spdep)
library(corrgram)
library(MuMIn)
library(visreg)
library(sjPlot)
library(car)
library(rsq)
library(caret)
library(relimp)
library(MASS)
library(spatialreg)
library(spind)
library(raster)
library(lattice)
library(nlme)
library(scales)
library(lme4)
library(piecewiseSEM)
library(performance)
library(lmtest)
library(dplyr)
library(plyr)
library(FactoMineR)
library(factoextra)
library(cluster)
library(ape)
library(vegan)
library(rdacca.hp)
library(ggord)
library(gridExtra)
library(grid)
library(ggpubr)
library(ggpattern)
library(ggalt)
library(ggrepel)
library(patchwork)
library(shades)
library(tidyverse)
### Load data
# ADNe presence/abscence
biodiv=read.csv("01_Analyses_teleo/00_data/matrice_teleo_totale.csv", row.names=1) %>%
t(.)

# Indicateurs
ind_ports =read.csv("indicators_ports_oct21_repassee.csv",row.names=1)

# milieu naturel
meta_nat <- read.csv("metadata_milieu_naturel.csv", header=T)
ind_nat <- read.csv("indicators_milieu_naturel.csv", header=T, row.names=1) %>%
  select(colnames(ind_ports))

# Combine the two datasets
ind <- rbind(ind_nat,ind_ports) %>%
  mutate(DP_B_ratio = log10(DP_B_ratio))

# metadata
meta=read.csv("metadata_tot.csv", sep=",",row.names=1) %>%
  column_to_rownames(var="code_spygen")

# merge indicateurs et metadata 
indmeta=merge(ind,meta,by=0) %>%
  filter(R>4) # enlever échantillons avec moins de 4 espèces
rownames(indmeta)=indmeta[,"Row.names"]

# combiner toutes les données
data_all=merge(biodiv,indmeta,by=0)
rownames(data_all)=data_all[,"Row.names"]
head(data_all)
dim(data_all)

## matrice des traits
traits <- read.table("Functional_data_corrected_20220124.csv", sep=",", header=T)

# vecteurs de noms d'espèces commerciales / non-commerciales / menacées
commercial <- traits %>%
  filter(all_commercial_level == 1) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  filter(Species %in% colnames(biodiv)) %>%
  pull(Species)

noncomm <- traits %>%
  filter(all_commercial_level == 0) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  filter(Species %in% colnames(biodiv)) %>%
  pull(Species)

threatened <- traits %>%
  filter(RedList >= 1) %>%
  mutate(Species = gsub(" ", "_", Species)) %>%
  filter(Species %in% colnames(biodiv)) %>%
  pull(Species)

#### faire la dbRDA
RDA_all=capscale(data_all[,c(2:186)] ~ port, data_all, dist="jaccard", add =TRUE)

#############################################################
## Plot RDA port/hors port + species
#############################################################
# get scores
site_scores <- scores(RDA_all)$sites     ## separating out the site scores
species_scores <- scores(RDA_all)$species %>% data.frame()   ## separating out the species
rsqr <- round(RsquareAdj(RDA_all)$adj.r.squared*100, 2) # r squared

# get most differentiated species along first axis (les espèces qui "tirent" la dbRDA = celles qui ont le plus d'importance pour différencier les sites)
quant50 <- quantile(abs(species_scores$CAP1), probs = c(0.50))
species_scores_diff50 <- species_scores[which(abs(species_scores$CAP1) > quant50["50%"]),]
# add colour variable by commercial/non-commercial/threatened
species_scores_diff50$col <- rep("other", nrow(species_scores_diff50))
species_scores_diff50$col[rownames(species_scores_diff50) %in% commercial]   <- "commercial"
species_scores_diff50$col[rownames(species_scores_diff50) %in% threatened] <- "threatened"
species_scores_diff50$col <- factor(species_scores_diff50$col, levels = c("threatened", "commercial", "other"))

# remplace les _ par espace dans nom d'espèces
rownames(species_scores) <- gsub("_", " ", rownames(species_scores))
rownames(species_scores_diff50) <- gsub("_", " ", rownames(species_scores_diff50))

# extract the percentage variability explained by axes
sumdbrda <- summary(RDA_all)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
MDS1 <- round(sumdbrda$cont$importance["Proportion Explained", "MDS1"]*100, 1)

# add metadata
site_scores_groups <- cbind(site_scores,select(data_all, port))

# vecteur de couleurs
mycol = c("#FCBBA1FF","#7FCDBBFF")

##########################
# plot in ggplot
#########################
### Represente les sites par catégories 
grda_sites <- ggplot(site_scores_groups, aes(x= CAP1, y = MDS1)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") + # ligne horizontqle
  geom_vline(xintercept = 0, lty = 2, col = "grey") + # ligne verticale
  geom_encircle(aes(group = port,linetype = port,fill= port), s_shape = 1, expand = 0,
                alpha = 0.6, show.legend = FALSE) + # hull area 
  geom_point(aes(pch = port, fill = port), cex = 3, col = "black") +
  ylim(-3.7, 2.4) + # limites axes Y
  scale_fill_manual(values = as.vector(mycol),
                    name = "", labels = c("Hors port", "Port")) +
  scale_shape_manual(values = c(22,21),
                    name = "", labels = c("Hors port", "Port")) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("MDS1 (", MDS1, "%)"),
       title = "") +
  theme_bw() +
  guides(fill=guide_legend(ncol=2)) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size =13, colour = "gray25"),
        axis.title = element_text(size = 14, colour = "gray25"),
        legend.position = c(1, 0),             # position in bottom right corner
        legend.justification = c(1, 0),        # correct legend justificaton
        legend.box.margin=margin(c(0,1,1,1)),  # add margin as to not overlap with axis box
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        legend.box.background=element_rect(fill='transparent', colour='transparent'),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) 
grda_sites

### Represente les espèces
grda_species <- ggplot() + 
  geom_segment(data= species_scores_diff50, aes(x=0, xend=CAP1,y = 0, yend=MDS1, col = col),
               arrow=arrow(length=unit(0.01,"npc")), size=0.5) + # most differentiated species
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  scale_color_manual(values = c("#e31a1c", "#1f78b4", "grey"), name = "", labels = c("Threatened species","Commercial species", "Other species")) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("MDS1 (", CAP2, "%)")) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 1.2))) +
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size =13, colour = "gray25"),
        axis.title = element_text(size = 14, colour = "gray25"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position = c(0, 1),              # position in top left corner
        legend.justification = c(0, 1),         # correct legend justificaton
        legend.box.margin=margin(c(2,2,2,2)),   # add margin as to not overlap with axis box
        legend.background = element_rect(fill =  alpha("white", 0.0)),
        legend.title = element_blank(),
        legend.text = element_text(size=18)) 
grda_species

# combine the 2 graphs
fig2 <-  (grda_sites + grda_species) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face="bold", size = 17))
fig2

# export figure
ggsave(plot = fig2, filename = "db-RDA_with_species.jpeg", 
       width = 15,  height = 7, dpi = 600)