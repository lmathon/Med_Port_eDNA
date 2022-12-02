library("VennDiagram")
library(ggvenn)
library(dplyr)
library(nVennR)
library(venn)
library(VennDiagram)
library(dplyr)
library(nVennR)
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


# load port data and extract species names
port <- read.csv("01_Analyses_teleo/00_data/teleo_presence.csv", row.names = 1)

port <- port[rowSums(port)!=0,]
port <- port[,colSums(port)!=0]

port <- rownames(port)
port <- gsub("_", " ", port)


# load outside data and extract species names
outside <- read.csv("01_Analyses_teleo/00_data/biodiv_milieu_naturel.csv", row.names = 1) 

outside <- outside[rowSums(outside)!=0,]
outside <- outside[,colSums(outside)!=0]

outside <- rownames(outside)
outside <- gsub("_", " ", outside)
  

# plot venn diagram simple
species_venn <- list(
  Ports = port,
  Hors_port = outside
)


vennplot <- ggvenn(species_venn, 
                         fill_color = c("#EE5E81", "#77CB81"),
                         stroke_size = 0.5, set_name_size = 4, show_percentage = TRUE, show_elements = FALSE)


# plot venn diagram proportional
species_venn2 <- createVennObj(nSets = 2, sNames = c('Hors_port', 'Ports'))
species_venn2 <- setVennRegion(species_venn2, c("Hors_port"), 32)
species_venn2 <- setVennRegion(species_venn2, c("Ports"), 35)
species_venn2 <- setVennRegion(species_venn2, c("Hors_port", "Ports"), 78)


vennplot2 <- plotVenn(nVennObj = species_venn2, systemShow=T)

####################################################################################
## Diagram with species lists
####################################################################################
## Venn diagram
venn.diagram(
  x = species_venn,
  category.names = c("Ports" , "Hors port"),
  cat.cex = 1.5, cex=1.5,
  col = c("#EE5E81", "#77CB81"),
  fill = c("#EE5E81", "#77CB81"),
  cat.col = c("#EE5E81", "#77CB81"),
  alpha = 0.2,
  filename="01_Analyses_teleo/03_Outputs/Species_Venn_Port_HorsPort1.png",
  output=F)

# Plot species names
png("01_Analyses_teleo/03_Outputs/Species_Venn_Port_HorsPort3.png", width = 1500, height = 1100)

# Setting a multi-panel graph 4x3
layout(
  matrix(1:5, ncol=5, byrow=TRUE),  # plot 12 graphs in 3 columns
  widths=c(1,1,1,1,1), # widths of each column
  heights=c(1)) # height of each row

# Print the names of the species that are only outside
plot(x = 0:1,                   # Create empty plot
     y = 0:1,
     ann = F,
     bty = "n",
     type = "n",
     xaxt = "n",
     yaxt = "n")
text(x = 0.55,                   # Add text to empty plot
     y = 0.5,
     paste(setdiff(outside, port)[1:50], collapse="\n"), 
     col= "black",
     cex = 2,
     font=3) # italic
plot(x = 0:1,                   # Create empty plot
     y = 0:1,
     ann = F,
     bty = "n",
     type = "n",
     xaxt = "n",
     yaxt = "n")
text(x = 0.55,                   # Add text to empty plot
     y = 0.5,
     paste(setdiff(outside, port)[51:100], collapse="\n"), 
     col= "black",
     cex = 2,
     font=3) # italic

# Print the names of the species that are only inside
plot(x = 0:1,                   # Create empty plot
     y = 0:1,
     ann = F,
     bty = "n",
     type = "n",
     xaxt = "n",
     yaxt = "n")
text(x = 0.55,                   # Add text to empty plot
     y = 0.5,
     paste(setdiff(port, outside), collapse="\n"),
     col = "black",
     cex = 2,
     font=3) # italic

dev.off()

#####################################################################################
## Venn diagram per category (reserve/lockdown)
#####################################################################################
# load port data and extract species names
adne <- read.csv("01_Analyses_teleo/00_data/teleo_presence_per_port.csv", row.names = 1)

port <- adne[rowSums(adne)!=0,]
port <- port[,colSums(port)!=0]

port <- rownames(port)
port <- gsub("_", " ", port)


# load outside data and extract species names
outside <- read.csv("01_Analyses_teleo/00_data/biodiv_milieu_naturel.csv", row.names = 1)
meta_nat <- read.csv2("00_Metadata/metadata_milieu_naturel.csv", header=T)
samples_reserve_N <- meta_nat %>% 
  filter(Confinement == "N" & protection == "reserve") %>%
  pull(SPYGEN_code)
samples_reserve_Y <- meta_nat %>% 
  filter(Confinement == "Y" & protection == "reserve") %>%
  pull(SPYGEN_code)
samples_out_N <- meta_nat %>% 
  filter(Confinement == "N" & protection == "outside") %>%
  pull(SPYGEN_code)
samples_out_Y <- meta_nat %>% 
  filter(Confinement == "Y" & protection == "outside") %>%
  pull(SPYGEN_code)

# Lockdown reserve
reserve_Y <- outside[,samples_reserve_Y]
reserve_Y <- reserve_Y[rowSums(reserve_Y)!=0,]
reserve_Y <- reserve_Y[,colSums(reserve_Y)!=0]

reserve_Y <- rownames(reserve_Y)
reserve_Y <- gsub("_", " ", reserve_Y)

# Not Lockdown reserve
reserve_N <- outside[,samples_reserve_N]
reserve_N <- reserve_N[rowSums(reserve_N)!=0,]
reserve_N <- reserve_N[,colSums(reserve_N)!=0]

reserve_N <- rownames(reserve_N)
reserve_N <- gsub("_", " ", reserve_N)

# Lockdown fished
fished_Y <- outside[,samples_out_Y]
fished_Y <- fished_Y[rowSums(fished_Y)!=0,]
fished_Y <- fished_Y[,colSums(fished_Y)!=0]

fished_Y <- rownames(fished_Y)
fished_Y <- gsub("_", " ", fished_Y)

# Not Lockdown fished
fished_N <- outside[,samples_out_N]
fished_N <- fished_N[rowSums(fished_N)!=0,]
fished_N <- fished_N[,colSums(fished_N)!=0]

fished_N <- rownames(fished_N)
fished_N <- gsub("_", " ", fished_N)


# plot venn diagram simple
species_venn <- list(
  Ports = port,
  Reserve_Unlock = reserve_N,
  Reserve_Lockdown = reserve_Y,
  Fished_Lockdown = fished_Y,
  Fished_Unlock = fished_N
)

## Venn diagram
venn.diagram(
  x = species_venn,
  height = 2500, 
  width = 3100,
  resolution = 300, imagetype = "tiff", 
  units = "px",
  # category.names = c("Ports" , "Hors port"),
  cat.cex = 1.2, cex=1,
  cat.just=list(c(0,2),c(0,-5),c(0,0),c(0,0),c(1,-4)),
  col = viridis(5),
  fill = viridis(5),
  cat.col = viridis(5),
  alpha = 0.2,
  filename="01_Analyses_teleo/03_Outputs/Species_Venn_categories.png",
  output=F)

#################################################################################
## Table species list
#################################################################################
reserve <- union(reserve_N,reserve_Y)
lockdown <- union(reserve_Y,fished_Y)
fished <- union(fished_Y, fished_N)
natural <- union(reserve, fished)
unlock <- union(reserve_N,fished_N)

port_only <- setdiff(port, natural)
reserve_only <- setdiff(reserve, union(fished,port))
lockdown_only <- setdiff(lockdown, union(unlock,port))

## Create table with species list
max_length <- max(c(length(port_only), length(reserve_only), length(lockdown_only)))    # Find out maximum length
max_length 

sp_table <- data.frame(Port_only = c(port_only,                 # Create data frame with unequal vectors
                            rep(NA, max_length - length(port_only))),
                   Reserve_only = c(reserve_only,
                            rep(NA, max_length - length(reserve_only))),
                   Lockdown_only = c(lockdown_only,
                                    rep(NA, max_length - length(lockdown_only))))
sp_table  
write.csv(sp_table, "01_Analyses_teleo/03_Outputs/Species_list_per_category.csv", row.names=F)

######### sp per port
sp_port <- adne[port_only,] %>%
  mutate(Nb_port = rowSums(., na.rm=TRUE))
write.csv(sp_port, "01_Analyses_teleo/03_Outputs/Species_list_per_port.csv")
