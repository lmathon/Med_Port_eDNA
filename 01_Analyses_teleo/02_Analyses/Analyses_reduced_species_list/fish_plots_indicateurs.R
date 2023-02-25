## Librairies
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(lemon)
library(magick)
library(mapdata)


## Load data
data <- read.csv("01_Analyses_teleo/00_data/reduced_matrice_teleo_port.csv", header=T, row.names=1) %>%
  filter(rowSums(.) > 0) %>%
  t(.)  %>%
  as.data.frame(.)

meta <- read.csv2("00_Metadata/metadata_port.csv", header=T)
biohut <- meta %>% 
  filter(type == "BIOHUT_port") %>%
  pull(code_spygen)

traits <- read.csv("01_Analyses_teleo/00_data/Functional_data_corrected_20220124.csv", header=T)

ind_ports <- read.csv("01_Analyses_teleo/00_data/Indicators_ports_2022_per_port_reduced_sp_list.csv", header=T, row.names=1) %>%
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


###############################################################################################
## Plot indicator values for ports
###############################################################################################

################################################################################
## Map of species richness for each port for the 2 campaigns
################################################################################
########################################################################################
## MAP
#########################################################################################
R_tot <- read.csv("01_Analyses_teleo/00_data/Richness_total_port_reduced_sp_list.csv") %>%
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

### Save the map
image_write(final, "01_Analyses_teleo/02_Analyses/Analyses_reduced_species_list/Outputs/Map_Richesse_total.png", density=300)

# Save in RData
save(final, file = "01_Analyses_teleo/04_Plots/WebFig4_a_Map_Richesse_total.RData")
