library(ggsn)
library(ggmap)
library(maps)
library(magick)
library(mapdata)
library(tidyselect)
library(scales)
library(dplyr)
library(radiant.data)
library(rsq)
library(purrr)



## Load data
data <- read.csv("02_Analyses_metazoa/00_data/metazoa_presence.csv")
data$scientific_name <- gsub("_", " ", data$scientific_name)

# make species names row names
rownames(data) <- data$scientific_name 
data <- data[,-1]


meta <- read.csv("00_Metadata/metadata_port.csv", header=T)

data2 <- data %>%
  t(.) %>%
  as.data.frame(.) %>%
  rownames_to_column(var="code_spygen") %>%
  left_join(meta, by="code_spygen") %>%
  # remove the biohut samples
  filter(habitat == "Port") %>%
  group_by(site, Campaign) %>%
  summarise_at(2:276, sum) %>% 
  # convert to presence/abscence
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  # create a new var combining site and campaign column
  mutate(Names = paste(site,Campaign, sep="_")) %>%
  column_to_rownames(var="Names") %>%
  select(3:ncol(.)) %>%
  t() %>%
  as.data.frame()

richness <- data.frame(site=colnames(data2),richness=apply(data2, 2, sum))


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
  labs(title = "Metazoa species richness") +
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
df <- richness %>%
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
  df_i <- df[which(df$Site %in% ports[i]),c("richness", "Site", "Campaign")]
  colnames(df_i) <- c("Y", "Site", "Sampling")
  to_plot <- df_i %>% 
    group_by(Sampling) %>% 
    summarise(value = mean(Y),
              sd=sd(Y)) 
  
  bp[[i]]<-ggplot(data=to_plot, aes(x=Sampling, y=value, fill=Sampling)) +
    geom_bar( position = 'dodge',stat="identity") +
    scale_fill_manual(values=c("lightsalmon", "lightblue2")) +
    ylim(0,110) +
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
image_write(final, "02_Analyses_metazoa/03_Outputs/Map_Richesse.png", density=300)
