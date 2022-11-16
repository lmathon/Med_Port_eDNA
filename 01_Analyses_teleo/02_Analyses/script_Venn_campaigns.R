library("VennDiagram")
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

## Load data
data <- read.csv("01_Analyses_teleo/00_data/teleo_presence.csv", header=T, row.names=1) %>%
  filter(rowSums(.) > 0) %>%
  t(.)  %>%
  as.data.frame(.) %>%
  rownames_to_column(var="code_spygen")

meta <- read.csv("00_Metadata/metadata_port.csv", header=T)

# list samples per campaign
Oct21_code <- meta %>%
  filter(Campaign == "October21") %>%
  filter(habitat == "Port") %>%
  filter(site != "Marseillan") %>%
  pull(code_spygen)

June22_code <- meta %>%
  filter(Campaign == "June22") %>%
  pull(code_spygen)

Oct21 <- data %>%
  filter(code_spygen %in% Oct21_code) %>%
  column_to_rownames(var="code_spygen") %>%
  t(.) %>%
  as.data.frame(.) %>%
  dplyr::filter(rowSums(.) > 0) %>%
  t(.)  %>%
  colnames(.)

June22 <- data %>%
  filter(code_spygen %in% June22_code) %>%
  column_to_rownames(var="code_spygen") %>%
  t(.) %>%
  as.data.frame(.) %>%
  dplyr::filter(rowSums(.) > 0) %>%
  t(.)  %>%
  colnames(.)

# plot venn diagram simple
species_venn <- list(
  Oct.21 = Oct21,
  Jun.22 = June22
)


## Venn diagram
venn.diagram(
  x = species_venn,
  category.names = c("Oct.21" , "Jun.22"),
  cat.cex = 1.5, cex=1.5,
  col = c("lightsalmon", "lightblue2"),
  fill = c("lightsalmon", "lightblue2"),
  cat.col = c("lightsalmon", "lightblue2"),
  alpha = 0.2,
  filename="01_Analyses_teleo/03_Outputs/Species_Venn_Campaign1.png",
  output=F)

# Plot species names
png("01_Analyses_teleo/03_Outputs/Species_Venn_Campaign3.png", width = 1200, height = 900)

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
     paste(setdiff(Oct21, June22), collapse="\n"), 
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
     paste(setdiff(June22, Oct21), collapse="\n"),
     col = "black",
     cex = 2,
     font=3) # italic

dev.off()
