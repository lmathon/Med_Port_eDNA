## Librairies
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(lemon)
library(magick)
library(mapdata)

# select species to keep
species_to_keep <- readRDS("01_Analyses_teleo/00_data/species_to_keep_ports.RDS")

## Load data
data <- read.csv("01_Analyses_teleo/00_data/teleo_presence.csv", header=T) %>%
  filter(scientific_name %in% species_to_keep) %>%
  column_to_rownames(var="scientific_name") %>%
  filter(rowSums(.) > 0) %>%
  t(.)  %>%
  as.data.frame(.)

meta <- read.csv2("00_Metadata/metadata_port.csv", header=T)
biohut <- meta %>% 
  filter(type == "BIOHUT_port") %>%
  pull(code_spygen)

traits <- read.csv("01_Analyses_teleo/00_data/Functional_data_corrected_20220124.csv", header=T)

ind_ports <- read.csv("01_Analyses_teleo/00_data/Indicators_ports_2022_per_filter.csv", header=T, row.names=1) %>%
  mutate(Location = "port") %>%
  dplyr::filter(!rownames(.) %in% biohut)

# Charger datas milieu naturel (= réserve et hors réserve)
meta_nat <- read.csv2("00_Metadata/metadata_milieu_naturel.csv", header=T)
ind_nat <- read.csv("01_Analyses_teleo/00_data/indicators_milieu_naturel.csv", header=T, row.names = 1) %>%
  rownames_to_column(var="Sample") %>%
  inner_join(meta_nat[,c("SPYGEN_code", "protection")], by=c("Sample"="SPYGEN_code")) %>%
  dplyr::rename(Location = protection) %>%
  column_to_rownames(var="Sample") %>%
  dplyr::select(colnames(ind_ports))

# Combiner les deux datasets
ind_all <- rbind(ind_nat,ind_ports) %>%
  mutate(DeBRa = log10(DeBRa))

# Boxplots
# Renommer les points milieu naturel = outside et reserve en "hors port"
ind_all$Location <- sub(pattern = "outside",  replacement = "hors port", ind_all$Location)
ind_all$Location <- sub(pattern = "reserve",  replacement = "hors port", ind_all$Location)

ind_all$Location <- factor(ind_all$Location,                                    # Factor levels in decreasing order
                            levels = c("hors port","port"))


p <-list()
## Draw the plot
l=1
for (i in c(1,4,8,11)) { # Choix des indicateurs à présenter
  # Gather data
  dat_i <- ind_all[,c(i,15)]
  colnames(dat_i)[1] <- "Y"
  
  p[[l]]<-ggplot(data=dat_i, aes(x=Location, y=Y, color=Location)) +
  geom_boxplot(notch=F) +
  scale_color_manual(values=c("darkblue", "cyan4",  "darkorange")) +
  ylab(colnames(ind_all)[i]) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.y=element_text(colour="black",size=16)) +
  theme(axis.text.x=element_text(colour="black",size=18)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_text(colour="black",size=18))+ 
  #stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size=4)
  stat_compare_means(aes(label = paste0("p = ", ..p.signif..)), size=6)

  l=l+1  
}

## Save plot
png("01_Analyses_teleo/03_Outputs/Boxplots_indicators_port_milieu_nat3.png", 
    width = 900, height = 1000) 
do.call(grid.arrange,c(p, list(ncol=2)))
dev.off()

################################################################################
####### Plot lockdown / reserve / ports
################################################################################
# load data
ind_ports <- read.csv("01_Analyses_teleo/00_data/Indicators_ports_2022_per_filter.csv", header=T, row.names=1) %>%
  dplyr::filter(!rownames(.) %in% biohut)

# Combiner les deux datasets
ind_all <- rbind(ind_nat[,1:14],ind_ports) %>%
  mutate(DeBRa = log10(DeBRa))

meta_tot <- read.csv("00_Metadata/metadata_tot.csv", header=T)

meta_tot <- meta_tot %>%
  filter(habitat != "BIOHUT_port")

ind_all <- ind_all %>%
  rownames_to_column(var="code_spygen") %>%
  inner_join(meta_tot, by="code_spygen") %>%
  mutate_at('Confinement', as.factor)

ind_all$habitat <- factor(ind_all$habitat, levels=c("reserve", "outside", "Port"))
ind_all$Confinement <- factor(ind_all$Confinement, labels=c("Unlock", "Lockdown"))

ind_names <- c("Total species richness", "Cryptobenthic species richness", "Threatened species richness", "Commercial species richness")  

p <-list()
## Draw the plot
l=1
for (i in c(1,4,8,11)) { # Choix des indicateurs à présenter
  # Gather data
  dat_i <- ind_all[,c(i,18,21)]
  colnames(dat_i)[1] <- "Y"
  
  p[[l]]<-ggplot(data=dat_i, aes(x=habitat, y=Y, fill=habitat)) +
    geom_boxplot(notch=F, alpha=0.5) +
    facet_wrap(~Confinement) +
    scale_fill_manual(values=c("darkblue", "cyan4",  "darkorange"),
                       labels=c("reserve", "fished", "port")) +
    ylab(ind_names[l]) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    theme(legend.position="none") +
    theme(axis.text.y=element_text(colour="black",size=16)) +
    theme(axis.text.x=element_text(colour="black",size=18)) +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_text(colour="black",size=18))+ 
    theme(strip.text = element_text(size = 18)) +
    #stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size=4)
    stat_compare_means(aes(label = paste0("p = ", ..p.signif..)), size=6)
  
  l=l+1  
}

## Save plot
png("01_Analyses_teleo/03_Outputs/Boxplots_indicators_per_category.png", 
    width = 1200, height = 1000) 
do.call(grid.arrange,c(p, list(ncol=2)))
dev.off()

####
## Tukey test between port and lockdown
ind_all <- ind_all %>%
  unite(group, c("habitat", "Confinement")) %>%
  mutate(group = as.factor(group))


# Create a matrix to store the coefficient and pvalues
compa <- c("outside_Unlock - outside_Lockdown" ,  "Port - outside_Lockdown"   ,   "reserve_Lockdown - outside_Lockdown",
            "reserve_Unlock - outside_Lockdown" ,  "Port - outside_Unlock" ,       "reserve_Lockdown - outside_Unlock",  
           "reserve_Unlock - outside_Unlock",     "reserve_Lockdown - Port" ,     "reserve_Unlock - Port",       
           "reserve_Unlock - reserve_Lockdown" )
res_tukey <- matrix(NA, 10, (2*ncol(Y)),
               dimnames= list(compa,
                              paste(rep(ind_names,each=2),c("coef", "pval"), sep="_" )))

# Make a loop to do the glm for each indicator
l=1

for (i in c(1,4,8,11)) { # Choix des indicateurs à présenter
  # Gather data
  dat_i <- ind_all[,c(i,18)]
  colnames(dat_i)[1] <- "Y"
  
  portAnova_i <- aov(Y ~ group, data = dat_i)
  portTukey_i <- TukeyHSD(portAnova_i, conf.level=.9) 
  
  res_tukey[,l:(l+1)] <- portTukey_i$group[,c(1,4)]
  
  l=l+2

}


write.csv(res_tukey,"01_Analyses_teleo/03_Outputs/PostHoc_test_indicators.csv")

###############################################################################################
## Plot indicator values for ports
###############################################################################################

################################################################################
## Map of species richness for each port for the 2 campaigns
################################################################################

########################################################################################
## MAP
#########################################################################################
ind_ports <- read.csv("01_Analyses_teleo/00_data/indicators_ports_2022_per_port.csv", header=T, row.names=1)
R_tot <- read.csv("01_Analyses_teleo/00_data/Richness_total_port.csv") %>%
  mutate(Campaign ="Total") %>%
  dplyr::rename(R = R_total,
         Site = Port) %>%
  dplyr::select(Site, Campaign, R)

# Download the map for the Mediterranean Sea
wH <- map_data("worldHires",  xlim=c(-8,37), ylim=c(29.5,47)) # subset polygons surrounding med sea

## create map and save
fig <- image_graph(width = 1200, height = 1100, res = 300)
ggplot() +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), 
               alpha=0.7,fill = "gray80", color = NA) +
  coord_fixed(xlim=c(3,7.3), ylim=c(42,44), ratio=1.2)+
  #geom_point(aes(x = longitude_start_DD, y = latitude_start_DD), data=df, pch=21, col = "black") +
  #scalebar(df, dist = 10, dist_unit = "km",
  #         transform = TRUE, model = "WGS84",
  #         st.bottom=T, location = "bottomright",
  #         anchor = c(x=8.95,y=42.495),
  #         st.dist = 0.05) +
  labs(title = "Species richness") +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        title = element_text(colour="black", size=7, face="bold"))
dev.off()


### put data together
df <- ind_ports %>%
  rownames_to_column(var="RowNames") %>%
  separate(RowNames, c("Site", "Campaign"), "_", extra = "merge") %>%
  mutate(Campaign = factor(Campaign, levels = c("October21", "June22"))) %>%
  dplyr::select(Site, Campaign, R) %>%
  bind_rows(R_tot) %>%
  mutate(Campaign = factor(Campaign, levels = c("Total","October21", "June22"))) 



# list of port names
ports <- unique(df$Site)
port_names <- ports
port_names[7] <- "St.Maries Mer"

## create list of the barplots
bp <- list()

for (i in 1: 7) {
  to_plot <- df[which(df$Site %in% ports[i]),]
  
  bp[[i]]<-ggplot(data=to_plot, aes(x=Campaign, y=R, fill=Campaign)) +
    geom_bar( position = 'dodge',stat="identity") +
    scale_fill_manual(values=c("lightblue","lightgreen", "sandybrown")) +
    ylim(0,78) +
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines 
      legend.position="none") +
    labs(y="Species richness",
         title=port_names[i]) +
    scale_x_discrete(labels=c("Total","Autumn", "Spring")) +
    theme(axis.text.y=element_text(colour="black",size=7)) +
    theme(axis.text.x=element_text(colour="black",size=7, angle=30, vjust=0.8)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x=element_blank())  +
    theme(plot.title = element_text(colour="black", size=7, face="bold")) 
  
}
names(bp) <- ports
bp

## create images with the barplots
barfig1 <- image_graph(width = 200, height = 220, res = 150, bg = 'transparent')
bp[[1]]
dev.off()

barfig2 <- image_graph(width = 200, height = 220, res = 150, bg = 'transparent')
bp[[2]]
dev.off()

barfig3 <- image_graph(width = 200, height = 220, res = 150, bg = 'transparent')
bp[[3]]
dev.off()

barfig4 <- image_graph(width = 200, height = 220, res = 150, bg = 'transparent')
bp[[4]]
dev.off()

barfig5 <- image_graph(width = 200, height = 220, res = 150, bg = 'transparent')
bp[[5]]
dev.off()

barfig6 <- image_graph(width = 210, height = 250, res = 150, bg = 'transparent')
bp[[6]]
dev.off()

barfig7 <- image_graph(width = 200, height = 220, res = 150, bg = 'transparent')
bp[[7]]
dev.off()

final <- image_composite(fig, barfig1, offset = "+100+460")
final <- image_composite(final, barfig2, offset = "+970+250") #
final <- image_composite(final, barfig3, offset = "+620+480") #
final <- image_composite(final, barfig4, offset = "+250+380")
final <- image_composite(final, barfig5, offset = "+820+540") #
final <- image_composite(final, barfig6, offset = "+100+650") #
final <- image_composite(final, barfig7, offset = "+440+330") #
final

# Save in RData
save(final, file = "01_Analyses_teleo/04_Plots/Fig2_a_Map_Richesse_total.RData")

### Save the map
image_write(final, "01_Analyses_teleo/03_Outputs/Map_Richesse_total.png", density=300)

#############################################################################################
## Barplot per port per campaign
#############################################################################################
### put data together
df <- ind_ports %>%
  rownames_to_column(var="RowNames") %>%
  separate(RowNames, c("Site", "Campaign"), "_", extra = "merge") %>%
  mutate(Campaign = factor(Campaign, levels = c("October21", "June22"))) %>%
  mutate(Site = factor(Site, levels = c("Cannes", "Porquerolles", "La Ciotat",
                                        "Saintes Maries de la Mer", 
                                        "Marseillan", "Agde", "Port Vendres")))


# list of port names
ports <- levels(df$Site)
port_names <- ports
port_names[4] <- "St.Maries Mer"

## vector of indicator names
names <- c("Species richness", "Functional Diversity", "Phylogenetic Diversity",
           "Cryptobenthic richness", "DeBRa", "Non-Indigenous species", "Threatened species",
           "Elasmobranch species", "Commercial species")

bp <- list()
e=1

for (i in c(3,4,15,6,7,9,10,11,13)) {
  to_plot <- df[,c(i,1,2)]
  colnames(to_plot) <- c("Y", "Site", "Sampling")
  
  
  bp[[e]]<-ggplot(data=to_plot, aes(x=Site, y=Y, fill=Sampling)) +
    geom_bar( position = 'dodge',stat="identity") +
    scale_fill_manual(values=c("lightsalmon", "lightblue2"),
                        name="Sampling",
                        labels = c("Oct. 21", "Jun. 22")) +
    labs(y=colnames(df)[i],
         title=names[e]) +
    scale_x_discrete(labels=port_names) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="bottom",
          axis.text.y=element_text(colour="black",size=8),
          axis.text.x=element_text(colour="black",size=9, angle=30, vjust=0.6),
          axis.title.y=element_text(colour="black",size=9),
          axis.title.x=element_blank(),
          title = element_text(colour="black", size=9, face="bold"),
          legend.text = element_text(colour="black", size=9),
          legend.title = element_text(colour="black", size=9,  face="bold")) 
  
  e=e+1

}

## Combine plots and save
Fig1 <- do.call("grid_arrange_shared_legend", c(bp, ncol=3, nrow=3, position="bottom"))
# save the figure
ggsave("01_Analyses_teleo/03_Outputs/Barplot_indicators_per_ports.jpeg",
       plot=Fig1,
       width = 31,
       height = 23,
       units = "cm")

#############################################################################################
## Boxplot per campaign
#############################################################################################
p <-list()
## Draw the plot
l=1
for (i in c(3,4,15,6,7,9,10,11,13)) { # Choix des indicateurs à présenter
  # Gather data
  dat_i <- df[,c(i,2)]
  colnames(dat_i)[1] <- "Y"
  
  p[[l]]<-ggplot(data=dat_i, aes(x=Campaign, y=Y, fill=Campaign)) +
    geom_boxplot(notch=F) +
    scale_fill_manual(values=c("lightsalmon", "lightblue2")) +
    labs(y=colnames(df)[i],
         title=names[l]) +
    scale_x_discrete(labels=c("Oct. 21","Jun. 22")) +
    theme_bw() +
    theme(legend.position="none") +
    theme(axis.text.y=element_text(colour="black",size=16)) +
    theme(axis.text.x=element_text(colour="black",size=18)) +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_text(colour="black",size=18))+ 
    #stat_compare_means(aes(label = paste0("p = ", ..p.format..)), size=4)
    stat_compare_means(aes(label = paste0("p = ", ..p.signif..)), size=6)
  
  l=l+1  
}

## Save plot
png("01_Analyses_teleo/03_Outputs/Boxplots_indicators_campaigns.png", 
    width = 1200, height = 1200) 
do.call(grid.arrange,c(p, list(ncol=3)))
dev.off()
