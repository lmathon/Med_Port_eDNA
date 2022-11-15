library(vegan)
library(ggplot2)
library(ade4)
library(magrittr) 
library(dplyr)  
library(dplyr)
library(rsq)
library(purrr)
library(ggplot2)
library(vegan)
library(betapart)

#TELE0 ESPECE

# load and prepare data : table with samples in rows and species in columns
metazoa <- read.csv("02_Analyses_metazoa/00_data/metazoa_reads.csv", sep=",", header = T, row.names = 1)

metazoa <- metazoa[,colSums(metazoa)!=0]

metazoa <- as.data.frame(t(metazoa))

# compute 'Bray' distance between samples
metazoa.bray <- vegdist(metazoa, method="bray")

# compute PCoA = ordination des points dans un espace
pcoa_metazoa <- dudi.pco(metazoa.bray, scannf = FALSE, nf = 3)

# select data to plot
metazoa2plot <- pcoa_metazoa$li

# ajout d'une colonne Ports et d'une colonne Habitat
meta <- read.csv("00_Metadata/metadata_port.csv", sep=",")

metazoa2plot$code_spygen <- rownames(metazoa2plot)
metazoa2plot <-left_join(metazoa2plot, meta[,c("code_spygen", "site", "habitat")])


##############################################################################################################################
## Load the eDNA data (matrix species per sample)
adne <- read.csv("02_Analyses_metazoa/00_data/metazoa_presence.csv", header=T, row.names = 1)

adne <- adne[rowSums(adne)!=0,]
adne <- adne[,colSums(adne)!=0]

sp_mat <- as.data.frame(t(adne))
sp_mat$Samples <- rownames(sp_mat)  
sp_mat <- inner_join(sp_mat, meta, by=c("Samples" = "code_spygen")) 

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


# plot

metazoa_plot <- ggplot(metazoa2plot, aes(x=A1, y=A2))+
  geom_point(aes(x=A1, y=A2, color=site, shape=habitat), size=2.5)+
  scale_shape_manual(values = c(16, 17), 
                     name = "Habitat",  labels = c("Biohut","Port"))+
  scale_color_manual(values = c("#FF575C","#F18F01","#FCF300", "#94F000", "#0087DB", "#A05CFF", "#F5A6E6"))+#, 
  #name = "Port", labels = c("Port_Vendres", "Agde", "Marseillan", "Saintes_Maries_de_la_Mer","La_Ciotat","Porquerolles","Cannes")) +
  scale_size_manual(values=c(2.5,2.5))+
  ylab("PCoA2")+
  xlab("PCoA1")+
  xlim(-0.6,0.8)+
  ylim(-0.6,0.4)+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "#FCFCFC", colour = NA),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill = NA),
        legend.text=element_text(size=15),
        axis.text=element_text(size=15),
        plot.title = element_text(size=18))+
  ggtitle("PCoA metazoa species composition among ports")

metazoa_plot

ggsave("02_Analyses_metazoa/03_Outputs/PCoA_metazoa_species.png", width = 11, height = 8)

