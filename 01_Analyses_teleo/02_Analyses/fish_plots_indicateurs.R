## Librairies
library(dplyr)
library(rsq)
library(fastDummies)
library(tidyr)
library(tidyverse)
library(purrr)
library(scales)
library(margins)
library(ggplot2)
library(gridExtra)
library(viridis)
library(ggpubr)
library(fishtree)


## Load data
data <- read.csv("01_Analyses_teleo/00_data/teleo_presence.csv", header=T, row.names=1) %>%
  filter(rowSums(.) > 0) %>%
  t(.)  %>%
  as.data.frame(.)

meta <- read.csv("00_Metadata/metadata_port.csv", header=T)

traits <- read.csv("01_Analyses_teleo/00_data/Functional_data_corrected_20220124.csv", header=T)

ind_ports <- read.csv("01_Analyses_teleo/00_data/indicators_ports_2022_per_port.csv", header=T, row.names=1) %>%
  mutate(Location = "port")

# Charger datas milieu naturel (= réserve et hors réserve)
meta_nat <- read.csv("00_Metadata/metadata_milieu_naturel.csv", header=T)
ind_nat <- read.csv("01_Analyses_teleo/00_data/indicators_milieu_naturel.csv", header=T, row.names=1) %>%
  inner_join(meta_nat[,c("SPYGEN_code", "protection")], by=c("Sample"="SPYGEN_code")) %>%
  rename(Location = protection) %>%
  column_to_rownames(var="Sample") %>%
  select(colnames(ind_ports))

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
for (i in c(1,2,8,11)) { # Choix des indicateurs à présenter
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

###############################################################################################
## Plot indicator values for ports
###############################################################################################

################################################################################
## Map of species richness for each port for the 2 campaigns
################################################################################
library(ggsn)
library(ggmap)
library(maps)
library(magick)
library(mapdata)
########################################################################################
## MAP
#########################################################################################

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
  mutate(Campaign = factor(Campaign, levels = c("October21", "June22")))

# list of port names
ports <- unique(df$Site)
port_names <- ports
port_names[7] <- "St.Maries Mer"

## create list of the barplots
bp <- list()

for (i in 1: 7) {
  df_i <- df[which(df$Site %in% ports[i]),c("R", "Site", "Campaign")]
  colnames(df_i) <- c("Y", "Site", "Sampling")
  to_plot <- df_i %>% 
    group_by(Sampling) %>% 
    summarise(value = mean(Y),
              sd=sd(Y)) 
  
  bp[[i]]<-ggplot(data=to_plot, aes(x=Sampling, y=value, fill=Sampling)) +
    geom_bar( position = 'dodge',stat="identity") +
    scale_fill_manual(values=c("lightsalmon", "lightblue2")) +
    ylim(0,56) +
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines 
      legend.position="none") +
    labs(y="Richness",
         title=port_names[i]) +
    scale_x_discrete(labels=c("Oct. 21","Jun. 22")) +
    theme(axis.text.y=element_text(colour="black",size=10)) +
    theme(axis.text.x=element_text(colour="black",size=10, angle=30, vjust=0.8)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x=element_blank())  +
    theme(plot.title = element_text(colour="black", size=9, face="bold", hjust=0.5)) 
  
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
final <- image_composite(final, barfig3, offset = "+600+440") #
final <- image_composite(final, barfig4, offset = "+250+380")
final <- image_composite(final, barfig5, offset = "+820+480") #
final <- image_composite(final, barfig6, offset = "+100+650") #
final <- image_composite(final, barfig7, offset = "+440+330") #
final

### Save the map
image_write(final, "01_Analyses_teleo/03_Outputs/Map_Richesse.png", density=300)
