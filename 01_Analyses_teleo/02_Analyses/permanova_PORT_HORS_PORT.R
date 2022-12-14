### Load the libraries
library(dplyr)
library(rsq)
library(purrr)
library(ggplot2)
library(vegan)
library(betapart)

#COMPARAISON PORT ET HORS PORT TOUT CONFONDU SANS LES BIOHUTS

# Metadata
meta_tot <- read.csv2("00_Metadata/metadata_tot.csv", header=T)

meta_tot <- meta_tot %>%
  filter(habitat != "BIOHUT_port")


## Load the eDNA data (matrix species per sample)
adne <- read.csv("01_Analyses_teleo/00_data/matrice_teleo_totale.csv", header=T, row.names=1)


# Enlever les pts biohuts
keep <- names(adne)[(names(adne) %in% meta_tot$code_spygen)]
adne <- adne[, keep]

adne <- adne[rowSums(adne)!=0,]
adne <- adne[,colSums(adne)!=0]


# Joindre matrice metadata et ADNe
sp_mat <- adne %>%
  t(.) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(var="Samples") %>%
  inner_join(meta_tot, by=c("Samples" = "code_spygen"))

## Create dissimilarity matrix (Jaccard distance for presence/absence)
sp.dist<-vegdist(t(adne), method='jaccard')

## Partition beta-diversity in turnover and nestedness
beta.dist <- beta.pair(t(adne), index.family = "jac")
turnover <- beta.dist[[1]]
nestedness <- beta.dist[[2]]

## Perform PERMANOVA : tester la significativite de la variable "habitat" (= port ou hors port)
# sur la diversite beta, le turnover et la nestedness
habitat.div<-adonis2(sp.dist~port, data=sp_mat, permutations = 999, method="jaccard")
habitat.div 
habitat.turn<-adonis2(turnover~port, data=sp_mat, permutations = 999, method="jaccard")
habitat.turn 
habitat.nest<-adonis2(nestedness~port, data=sp_mat, permutations = 999, method="jaccard")
habitat.nest 

## Multivariate dispersion : tester l'homogénéité des variances entre les habitats
dispersion_m<-betadisper(sp.dist, group=sp_mat$port)
permutest(dispersion_m)
dispersion_turn<-betadisper(turnover, group=sp_mat$port)
permutest(dispersion_turn)
dispersion_nest<-betadisper(nestedness, group=sp_mat$port)
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
grp1.hull<-seg.data[seg.data$group=="true",1:3][chull(seg.data[seg.data$group=="true",2:3]),]
grp2.hull<-seg.data[seg.data$group=="false",1:3][chull(seg.data[seg.data$group=="false",2:3]),]
all.hull<-rbind(grp1.hull,grp2.hull)

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
  legend.title = element_text(size = 18),
  legend.text = element_text(size = 18),
  legend.key = element_blank())

# plot
panel.a<-ggplot() + 
  geom_polygon(data=all.hull,aes(x=v.PCoA1,y=v.PCoA2, colour=group, fill=after_scale(alpha(colour, 0.3)))) +
  geom_point(data=centroids[,1:3], aes(x=PCoA1,y=PCoA2,colour=grps, shape=grps),size=4) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, colour=group, shape=group),size=2) +
  theme_mine + 
  scale_colour_manual(values=c("#FCBBA1FF","#7FCDBBFF"), labels = c("Sea","Port")) +
  scale_shape_manual(values=c(19,17), labels = c("Sea","Port")) +
  labs(title="",x="PCoA 1",y="PCoA 2",
       colour = "Sampling", shape = "Sampling")
panel.a


ggsave("01_Analyses_teleo/03_Outputs/PCoA_horsport_port.png", width = 11, height = 8)




