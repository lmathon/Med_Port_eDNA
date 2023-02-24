library(betapart)
library(tidyverse)
library(vegan)
library(ggalt)

### Load data
# eDNA presence/abscence
biodiv=read.csv("01_Analyses_teleo/00_data/matrice_teleo_totale.csv", row.names=1) %>%
  t(.) %>%
  as.data.frame(.)

biodiv$code_spygen <- rownames(biodiv)

# metadata
meta=read.csv("00_Metadata/metadata_tot.csv", sep=";") 
meta <- meta %>%
  filter(code_spygen %in% biodiv$code_spygen)

# combiner data and metadata
data_all=left_join(biodiv,meta, by = 'code_spygen', multiple = 'first')
data_all <- data_all %>%
  distinct(code_spygen, .keep_all=T)
rownames(data_all)=data_all[,"code_spygen"]

data_all <- data_all %>%
  filter(!is.na(port))

# filter out lockdown (facultatif)
data_all <- data_all %>%
  filter(Confinement=="N")

# data for dbRDA
data_dbrda <- data_all[,c(1:289)]
meta_dbrda <- data_all[,c(290:ncol(data_all))]


# Compute beta_diversity
b <- betapart.core(data_dbrda)
beta <- beta.pair(b, "jaccard")
beta_tot <- beta$beta.jac
beta_tur <- beta$beta.jtu
beta_nes <- beta$beta.jne


########################################################################################################################
#### FUNCTION PLOT GRDA
plot_grda_sites <- function(df, fill_values = c("#7FC5C5","#FFC47E", "#7F7FC5")){
  grda_sites <- ggplot(df, aes(x= CAP1, y = CAP2)) +
    geom_hline(yintercept = 0, lty = 2, col = "grey") + # ligne horizontqle
    geom_vline(xintercept = 0, lty = 2, col = "grey") + # ligne verticale
    geom_encircle(aes(group = habitat,linetype = habitat,fill= habitat), s_shape = 1, expand = 0,
                  alpha = 0.6, show.legend = FALSE) + # hull area 
    geom_point(aes(pch = habitat, fill = habitat), cex = 3, col = "black") +
    #ylim(-3.7, 2.4) + # limites axes Y
    scale_fill_manual(values = as.vector(fill_values),
                      name = "", labels = c("Fished", "Port", "Reserve")) +
    scale_shape_manual(values = c(22,21,23),
                       name = "", labels = c("Fished", "Port", "Reserve")) +
    labs(x = paste0("CAP1 (", CAP1, "%)"), y = paste0("CAP2 (", CAP2, "%)"),
         title = "") +
    theme_bw() +
    guides(fill=guide_legend(ncol=3)) +
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
}

#### Compute dbRDA on total beta-diversity ####
RDA_all<- capscale(beta_tot ~ habitat + Condition(Longitude), meta_dbrda)

### Test the significance
anova(RDA_all, perm=999)
anova(RDA_all, perm=999, by="margin") # significant predictors : 
anova(RDA_all, perm=999, by="axis")

# get scores
site_scores <- scores(RDA_all)$sites     ## separating out the site scores
species_scores <- scores(RDA_all)$species %>% data.frame()   ## separating out the species
rsqr <- round(RsquareAdj(RDA_all)$adj.r.squared*100, 2) # r squared

# extract the percentage variability explained by axes
sumdbrda <- summary(RDA_all)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(meta_dbrda$code_spygen), rownames(site_scores)) # verify that data in same order
site_scores_groups <- cbind(site_scores,select(meta_dbrda, habitat))

### plot sites 
grda_sites <- plot_grda_sites(df = site_scores_groups)
grda_sites

# export figure
ggsave(plot = grda_sites, filename = "01_Analyses_teleo/03_Outputs/dbRDA_totale.jpeg", 
       width = 7,  height = 7, dpi = 600)

########################################################################################################################
#### Compute dbRDA on turnover beta-diversity ####
RDA_tur<- capscale(beta_tur ~ habitat + Condition(Longitude), meta_dbrda)

# get scores
site_scores <- scores(RDA_tur)$sites     ## separating out the site scores
species_scores <- scores(RDA_tur)$species %>% data.frame()   ## separating out the species
rsqr <- round(RsquareAdj(RDA_tur)$adj.r.squared*100, 2) # r squared

# extract the percentage variability explained by axes
sumdbrda <- summary(RDA_tur)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(meta_dbrda$code_spygen), rownames(site_scores)) # verify that data in same order
site_scores_groups <- cbind(site_scores,select(meta_dbrda, habitat))

### plot sites
grda_sites <- plot_grda_sites(df = site_scores_groups)
grda_sites

# export figure
ggsave(plot = grda_sites, filename = "01_Analyses_teleo/03_Outputs/dbRDA_totale_turnover.jpeg", 
       width = 7,  height = 7, dpi = 600)


########################################################################################################################
#### Compute dbRDA on nestedness beta-diversity ####
RDA_nes<- capscale(beta_nes ~ habitat + Condition(Longitude), meta_dbrda)

# get scores
site_scores <- scores(RDA_nes)$sites     ## separating out the site scores
species_scores <- scores(RDA_nes)$species %>% data.frame()   ## separating out the species
rsqr <- round(RsquareAdj(RDA_nes)$adj.r.squared*100, 2) # r squared

# extract the percentage variability explained by axes
sumdbrda <- summary(RDA_nes)
CAP1 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP1"]*100, 1)
CAP2 <- round(sumdbrda$cont$importance["Proportion Explained", "CAP2"]*100, 1)

# add metadata
identical(as.character(meta_dbrda$code_spygen), rownames(site_scores)) # verify that data in same order
site_scores_groups <- cbind(site_scores,select(meta_dbrda, habitat))

### plot sites 
grda_sites <- plot_grda_sites(df = site_scores_groups)
grda_sites

# export figure
ggsave(plot = grda_sites, filename = "01_Analyses_teleo/03_Outputs/dbRDA_totale_nestedness.jpeg", 
       width = 7,  height = 7, dpi = 600)