### Load the libraries
library(dplyr)
library(rsq)
library(purrr)
library(ggplot2)
library(vegan)
library(betapart)

#COMPARAISON PORT ET HORS PORT TOUT CONFONDU SANS LES BH

# Metadata
meta_tot <- read.csv("00_Metadata/metadata_tot.csv", header=T)

## Load the eDNA data (matrix species per sample)
adne <- read.csv("03_Analyses_eukaryotes/00_data/matrice_eukaryote_totale.csv", header=T, row.names=1)

adne <- adne[rowSums(adne)!=0,]
adne <- adne[,colSums(adne)!=0]

# Joindre matrice metadata et ADNe
sp_mat <- adne %>%
  t(.) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var="Samples") %>%
  inner_join(meta_tot, by=c("Samples" = "code_spygen"))

## create dissimilarity matrix (Jaccard distance for presence/absence)
sp.dist<-vegdist(t(adne), method='jaccard')

## Partition beta-diversity in turnover and nestedness
beta.dist <- beta.pair(t(adne), index.family = "jac")
turnover <- beta.dist[[1]]
nestedness <- beta.dist[[2]]

## perform PERMANOVA : tester la significativite de la variable "habitat"
# sur la diversite beta, le turnover et la nestedness
habitat.div<-adonis2(sp.dist~habitat, data=sp_mat, permutations = 999, method="jaccard")
habitat.div # signif
habitat.turn<-adonis2(turnover~habitat, data=sp_mat, permutations = 999, method="jaccard")
habitat.turn # signif
habitat.nest<-adonis2(nestedness~habitat, data=sp_mat, permutations = 999, method="jaccard")
habitat.nest # not signif

## Multivariate dispersion : tester l'homogeneite des variances entre les habitats
dispersion_m<-betadisper(sp.dist, group=sp_mat$habitat)
permutest(dispersion_m)
dispersion_turn<-betadisper(turnover, group=sp_mat$habitat)
permutest(dispersion_turn)
dispersion_nest<-betadisper(nestedness, group=sp_mat$habitat)
permutest(dispersion_nest)


###########################################################
## Plots with ggplot
###########################################################
library(gridExtra)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(grps=rownames(dispersion_m$centroids),data.frame(dispersion_m$centroids))
vectors<-data.frame(group=dispersion_m$group,data.frame(dispersion_m$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")


# create the convex hulls of the outermost points
grp1.hull<-seg.data[seg.data$group=="outside",1:3][chull(seg.data[seg.data$group=="outside",2:3]),]
grp2.hull<-seg.data[seg.data$group=="reserve",1:3][chull(seg.data[seg.data$group=="reserve",2:3]),]
grp3.hull<-seg.data[seg.data$group=="Port",1:3][chull(seg.data[seg.data$group=="Port",2:3]),]
grp4.hull<-seg.data[seg.data$group=="BIOHUT_port",1:3][chull(seg.data[seg.data$group=="BIOHUT_port",2:3]),]
all.hull<-rbind(grp1.hull,grp2.hull,grp3.hull,grp4.hull)

## Draw plot
# Create a theme
theme_mine <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 15),
  axis.text = element_text(size = 16, colour = "gray25"),
  axis.title = element_text(size = 18, colour = "gray25"),
  legend.position = "right",
  legend.title = element_blank(),
  legend.text = element_text(size = 18),
  legend.key = element_blank())

# plot
panel.a<-ggplot() + 
  geom_polygon(data=all.hull,aes(x=v.PCoA1,y=v.PCoA2, colour=group, fill=after_scale(alpha(colour, 0.3)))) +
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2,colour=grps),size=4) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, colour=group),size=2) +
  theme_mine + 
  scale_colour_manual(values=c("#F5A6E6","#F18F01","#0087DB"), labels = c("Biohut","Port", "Reserve")) +
  labs(title="PCoA teleost species among habitats",x="PCoA 1",y="PCoA 2")

panel.a


ggsave("03_Analyses_eukaryotes/03_Outputs/PCoA_habitats.png", width = 11, height = 8)




