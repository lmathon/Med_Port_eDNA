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

port <- as.data.frame(rownames(port))
names(port) <- "species"
  
# keep only species genus_species
port <- port %>%
  filter(grepl('_', species)) %>%
  filter(!grepl('_sp.', species)) 

port <- as.vector(port$species)


# load outside data and extract species names
outside <- read.csv("01_Analyses_teleo/00_data/biodiv_milieu_naturel.csv", row.names = 1) %>%
  t(.) %>%
  as.data.frame(.)

outside <- outside[rowSums(outside)!=0,]
outside <- outside[,colSums(outside)!=0]

outside <- rownames(outside)
  


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
png("01_Analyses_teleo/03_Outputs/Species_Venn_Port_HorsPort3.png", width = 1200, height = 900)

# Setting a multi-panel graph 4x3
layout(
  matrix(1:4, ncol=4, byrow=TRUE),  # plot 12 graphs in 3 columns
  widths=c(1,1,1,1), # widths of each column
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
     paste(setdiff(outside, port), collapse="\n"), 
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
