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

### Load data
# eDNA presence/abscence
biodiv=read.csv("01_Analyses_teleo/00_data/matrice_teleo_totale.csv", row.names=1) %>%
t(.) %>%
as.data.frame(.)

biodiv$code_spygen <- rownames(biodiv)


# Indicators
ind_ports =read.csv("01_Analyses_teleo/00_data/Indicators_ports_2022_per_filter.csv",row.names=1)

# milieu naturel
meta_nat <- read.csv2("00_Metadata/metadata_milieu_naturel.csv", header=T)
ind_nat <- read.csv("01_Analyses_teleo/00_data/indicators_milieu_naturel.csv", header=T, row.names = 1) %>%
  select(colnames(ind_ports))

# Combine the two datasets
ind <- rbind(ind_nat,ind_ports) 
ind$code_spygen <- rownames(ind)

# metadata
meta=read.csv("00_Metadata/metadata_tot.csv", sep=";") 

# merge indicateurs et metadata
indmeta <- left_join(ind, meta)
indmeta <- indmeta %>%
  distinct(code_spygen, .keep_all=T)
rownames(indmeta)=indmeta[,"code_spygen"]

# combiner toutes les données
data_all=left_join(biodiv,indmeta)
rownames(data_all)=data_all[,"code_spygen"]
head(data_all)
dim(data_all)

data_all <- data_all %>%
  filter(!is.na(port))

## matrice des traits
traits <- read.table("01_Analyses_teleo/00_data/Functional_data_corrected_20220124.csv", sep=",", header=T)

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
data_dbrda <- data_all[,c(1:184)]
meta_dbrda <- data_all[,c(185:ncol(data_all))]

# Compute beta_diversity
b <- betapart.core(data_dbrda)
beta <- beta.pair(b, "jaccard")
beta_tot <- beta$beta.jac
beta_tur <- beta$beta.jtu
beta_nes <- beta$beta.jne



RDA_all<- capscale(data_dbrda ~ port + Condition(Longitude), meta_dbrda, distance ="jaccard")

# get scores
site_scores <- scores(RDA_all)$sites     ## separating out the site scores
species_scores <- scores(RDA_all)$species %>% data.frame()   ## separating out the species
rsqr <- round(RsquareAdj(RDA_all)$adj.r.squared*100, 2) # r squared

# get most differentiated species along first axis (les espèces qui "tirent" la dbRDA = celles qui ont le plus d'importance pour différencier les sites)
quant50_cap1 <- quantile(abs(species_scores$CAP1), probs = c(0.5))
quant50_cap2 <- quantile(abs(species_scores$MDS1), probs = c(0.5))
quant50 <- rbind(quant50_cap1, quant50_cap2 )
species_scores_diff50_cap1 <- species_scores[which(abs(species_scores$CAP1) > quant50_cap1["50%"]),]
species_scores_diff50_cap2 <- species_scores[which(abs(species_scores$MDS1) > quant50_cap2["50%"]),]
species_scores_diff50 <- rbind(species_scores_diff50_cap1, species_scores_diff50_cap2)
species_scores_diff50 <- unique(species_scores_diff50)


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
identical(as.character(meta_dbrda$code_spygen), rownames(site_scores)) # verify that data in same order
site_scores_groups <- cbind(site_scores,select(meta_dbrda, port))

# vecteur de couleurs
mycol = c("#7FC5C5","#FFC47E")

# Plot sites 
grda_sites <- ggplot(site_scores_groups, aes(x= CAP1, y = MDS1)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey") + # ligne horizontqle
  geom_vline(xintercept = 0, lty = 2, col = "grey") + # ligne verticale
  geom_encircle(aes(group = port,linetype = port,fill= port), s_shape = 1, expand = 0,
                alpha = 0.6, show.legend = FALSE) + # hull area 
  geom_point(aes(pch = port, fill = port), cex = 3, col = "black") +
  ylim(-3.7, 2.4) + # limites axes Y
  scale_fill_manual(values = as.vector(mycol),
                    name = "", labels = c("Outside", "Port")) +
  scale_shape_manual(values = c(22,21),
                    name = "", labels = c("Outside", "Port")) +
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

# plot species
grda_species <- ggplot() + 
  geom_segment(data= species_scores_diff50, aes(x=0, xend=CAP1,y = 0, yend=MDS1, col = col),
               arrow=arrow(length=unit(0.01,"npc")), size=0.5) + # most differentiated species
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  scale_color_manual(values = c("#e31a1c", "#1f78b4", "grey"), name = "", labels = c("Threatened species","Commercial species", "Other species")) +
  labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("MDS1 (", MDS1, "%)")) +
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
fig2 <-  ggarrange(grda_sites,grda_species, ncol=2) 
fig2

# export figure
ggsave(plot = fig2, filename = "01_Analyses_teleo/03_Outputs/dbRDA_with_species.jpeg", 
       width = 15,  height = 7, dpi = 600)



