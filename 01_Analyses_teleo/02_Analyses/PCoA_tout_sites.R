### Load the libraries
library(dplyr)
library(rsq)
library(purrr)
library(ggplot2)
library(vegan)
library(betapart)
library(ade4)

#COMPARAISON PORT ET HORS PORT TOUT CONFONDU SANS LES BH

# Metadata
meta_tot <- read.csv("00_Metadata/metadata_tot.csv", header=T)

## Load the eDNA data (matrix species per sample)
adne <- read.csv("01_Analyses_teleo/00_data/matrice_teleo_totale.csv", header=T, row.names=1)


teleo <- as.data.frame(t(adne))

# compute 'Bray' distance between samples
teleo.bray <- vegdist(teleo, method="bray")

# compute PCoA = ordination des points dans un espace
pcoa_teleo <- dudi.pco(teleo.bray, scannf = FALSE, nf = 3)

# select data to plot
teleo2plot <- pcoa_teleo$li

# ajout d'une colonne Ports et d'une colonne Habitat
teleo2plot$code_spygen <- rownames(teleo2plot)
teleo2plot <-left_join(teleo2plot, meta_tot[,c("code_spygen", "site", "habitat")])

################################################################################################################################???
adne <- adne[rowSums(adne)!=0,]
adne <- adne[,colSums(adne)!=0]

# Joindre matrice metadata et ADNe
sp_mat <- adne %>%
  t(.) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var="Samples") %>%
  left_join(meta_tot, by=c("Samples" = "code_spygen"))

## create dissimilarity matrix (Jaccard distance for presence/absence)
sp.dist<-vegdist(t(adne), method='jaccard')

## Partition beta-diversity in turnover and nestedness
beta.dist <- beta.pair(t(adne), index.family = "jac")
turnover <- beta.dist[[1]]
nestedness <- beta.dist[[2]]



## perform PERMANOVA : tester la significativite de la variable "site" 
# sur la diversite beta, le turnover et la nestedness
habitat.div<-adonis2(sp.dist~site, data=sp_mat, permutations = 999, method="jaccard")
habitat.div 
habitat.turn<-adonis2(turnover~site, data=sp_mat, permutations = 999, method="jaccard")
habitat.turn 
habitat.nest<-adonis2(nestedness~site, data=sp_mat, permutations = 999, method="jaccard")
habitat.nest 

## Multivariate dispersion : tester l'homogeneite des variances entre les site
dispersion_m<-betadisper(sp.dist, group=sp_mat$site)
permutest(dispersion_m)
dispersion_turn<-betadisper(turnover, group=sp_mat$site)
permutest(dispersion_turn)
dispersion_nest<-betadisper(nestedness, group=sp_mat$site)
permutest(dispersion_nest)



##############################################################################################################################
library(scales)                                     # Amount of default colors
hex_codes1 <- hue_pal()(14)                             # Identify hex codes
hex_codes1

colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
            "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "grey80", "black",
            "#3F3C05", "#1D6F61", "#FE065E", "#FB61D7", "#FF66A8")
# plot

teleo_plot <- ggplot(teleo2plot, aes(x=A1, y=A2))+
  geom_point(aes(x=A1, y=A2, color=site), size=2.5)+
  scale_color_manual(values=colors)+#, 
  #name = "Port", labels = c("Port_Vendres", "Agde", "Marseillan", "Saintes_Maries_de_la_Mer","La_Ciotat","Porquerolles","Cannes")) +
  scale_size_manual(values=c(2.5,2.5))+
  ylab("PCoA2")+
  xlab("PCoA1")+
  xlim(-0.6,0.8)+
  ylim(-0.7,0.5)+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "#FCFCFC", colour = NA),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill = NA),
        legend.text=element_text(size=15),
        axis.text=element_text(size=15),
        plot.title = element_text(size=18))+
  ggtitle("PCoA teleost species composition among ports")

teleo_plot

ggsave("01_Analyses_teleo/03_Outputs/PCoA_teleo_tous_sites.png", width = 11, height = 8)
