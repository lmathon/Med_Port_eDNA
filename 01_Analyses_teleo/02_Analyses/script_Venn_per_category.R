## ---------------------------
##
## Venn Diagram:
##   1. Port, Reserve, Fished
##   2. Port, Lockdown, Unlock
##
## Date Created: 2023-02-07
## ----------------------------

#_libs
library(dplyr)
library(VennDiagram)
library(cowplot)

# Input
## Fun
rm_null_rowcol <- function(data){
    data <- data[rowSums(data) != 0,]
    data[,colSums(data) != 0]
}
## Port data
port_adne <- read.csv("01_Analyses_teleo/00_data/matrice_teleo_port.csv", row.names = 1)
port_adne <- rm_null_rowcol(port_adne)

## Outside data (Reserve + Fished)
outside <- read.csv("01_Analyses_teleo/00_data/biodiv_milieu_naturel.csv", row.names = 1)
meta_nat <- read.csv2("00_Metadata/metadata_milieu_naturel.csv", header = TRUE)

## Samples Reserve
### Unlock
samples_reserve_N <- meta_nat %>% 
  dplyr::filter(Confinement == "N" & protection == "reserve") %>%
  dplyr::pull(SPYGEN_code)
### Lockdown
samples_reserve_Y <- meta_nat %>% 
  dplyr::filter(Confinement == "Y" & protection == "reserve") %>%
  dplyr::pull(SPYGEN_code)
## Samples Fished
### Unlock
samples_out_N <- meta_nat %>% 
  dplyr::filter(Confinement == "N" & protection == "outside") %>%
  dplyr::pull(SPYGEN_code)
### Lockdown
samples_out_Y <- meta_nat %>% 
  dplyr::filter(Confinement == "Y" & protection == "outside") %>%
  dplyr::pull(SPYGEN_code)

## Reserve data
### Unlock
reserve_N <- outside[,samples_reserve_N]
reserve_N <- rm_null_rowcol(reserve_N)
### Lockdown
reserve_Y <- outside[ , samples_reserve_Y]
reserve_Y <- rm_null_rowcol(reserve_Y)

## Fished Data
### Unlock
fished_N <- outside[,samples_out_N]
fished_N <- rm_null_rowcol(fished_N)
### Lockdown
fished_Y <- outside[,samples_out_Y]
fished_Y <- rm_null_rowcol(fished_Y)

## Species list
### Port
port <- rownames(port_adne)
port <- gsub("_", " ", port)
### Reserve
reserve_N <- rownames(reserve_N)
reserve_N <- gsub("_", " ", reserve_N)
reserve_Y <- rownames(reserve_Y)
reserve_Y <- gsub("_", " ", reserve_Y)
### Fished
fished_N <- rownames(fished_N)
fished_N <- gsub("_", " ", fished_N)
fished_Y <- rownames(fished_Y)
fished_Y <- gsub("_", " ", fished_Y)

# Venn Diagram
## Fun
plot_venn <- function(list_cate, file_name, pal_fill, pal_col, cat_just){
  VennDiagram::venn.diagram(
    x = list_cate,
    filename = file_name,
    output=TRUE,
    
    # Set Output
    imagetype = "tiff" ,
    height = 2500 , 
    width = 3100 , 
    resolution = 300,
    unity = "px",
    
    # Set Circles
    fill = pal_fill,
    col = pal_col,
    lwd = 3,
    
    # Set Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set Categories
    cat.cex = 1.2,
    cat.just = cat_just,
    cat.fontfamily = "sans",
    cat.fontface = "bold",
    cat.col = pal_col
  )
}

## Venn diagram 1 - CATEGORY = Port, Reserve Unlock, Reserve Lockdown, Fished Unlock, Fished Lockdown
species_venn <- list(
  Port = port,
  Reserve_Unlock = reserve_N,
  Reserve_Lockdown = reserve_Y,
  Fished_Lockdown = fished_Y,
  Fished_Unlock = fished_N
)
plot_venn(list_cate = species_venn, file_name = "01_Analyses_teleo/03_Outputs/Species_Venn_categories.tiff", 
          pal_fill = c("#fdd8a8", "#b3b1d9", "#b3b1d9", "#b2dadb", "#b2dadb"), 
          pal_col = c("#fd9435", "#1c0087", "#b3b1d9", "#b2dadb", "#138889"),
          cat_just = list(c(.5, 2), c(0, -5), c(0, 0), c(.3, 0), c(1, -4)))

## Venn diagram 2 - CATEGORY = Port, Lockdown, Unlock
species_venn2 <- list(
  Port = port,
  Unlock = union(reserve_N, fished_N),
  Lockdown = union(reserve_Y, fished_Y)
)
plot_venn(list_cate = species_venn2, file_name = "01_Analyses_teleo/03_Outputs/Species_Venn_port_lockdown_unlock.tiff", 
          pal_fill = c("#fdd8a8", "#B2C6DA", "#B2C6DA"), 
          pal_col = c("#fd9435", "#6763b0", "#B2C6DA"),
          cat_just = list(c(0, 0), c(1, 0), c(.4, 0)))

## Venn diagram 3 - CATEGORY = Port, Reserve, Fished
species_venn3 <- list(
  Port = port,
  Reserve = union(reserve_N, reserve_Y),
  Fished = union(fished_N, fished_Y)
)
plot_venn(list_cate = species_venn3, file_name = "01_Analyses_teleo/03_Outputs/Species_Venn_port_reserve_fished.tiff", 
          pal_fill = c("#fdd8a8", "#b3b1d9", "#b2dadb"), 
          pal_col = c("#fd9435", "#1c0087", "#138889"),
          cat_just = list(c(0, 0), c(1, 0), c(.4, 0)))

## Save plots in publication

venn_portlockun <- plot_venn(list_cate = species_venn2, file_name = NULL, 
                pal_fill = c("#fdd8a8", "#B2C6DA", "#B2C6DA"), 
                pal_col = c("#fd9435", "#6763b0", "#B2C6DA"),
                cat_just = list(c(0, 0), c(.7, 0), c(.4, .2)))
venn_portreservfish <- plot_venn(list_cate = species_venn3, file_name = NULL, 
                pal_fill = c("#fdd8a8", "#b3b1d9", "#b2dadb"), 
                pal_col = c("#fd9435", "#1c0087", "#138889"),
                cat_just = list(c(0, 0), c(.7, 0), c(.4, .2)))

save(venn_portlockun, venn_portreservfish, file = "01_Analyses_teleo/04_Plots/Fig3_ab_venn.RData")

# Output table
reserve <- union(reserve_N,reserve_Y)
lockdown <- union(reserve_Y,fished_Y)
fished <- union(fished_Y, fished_N)
natural <- union(reserve, fished)
unlock <- union(reserve_N,fished_N)
## Only table
port_only <- setdiff(port, natural)
reserve_only <- setdiff(reserve, union(fished,port))
lockdown_only <- setdiff(lockdown, union(unlock,port))
only_list <- list(Port_only = port_only, Reserve_only = reserve_only, Lockdown_only = lockdown_only)
max_length <- max(sapply(only_list, length))
sp_table <- sapply(only_list, function(x){
    c(x, rep(NA, max_length - length(x)))
}) %>% data.frame()
write.csv(sp_table, "01_Analyses_teleo/03_Outputs/Species_list_per_category.csv", row.names=F)
## sp per port
sp_port <- port_adne[port_only,] %>%
  dplyr::mutate(Nb_port = rowSums(., na.rm=TRUE))
write.csv(sp_port, "01_Analyses_teleo/03_Outputs/Species_list_per_port.csv")
## sp in reserve
write.csv(list('Species' = reserve), "01_Analyses_teleo/03_Outputs/Species_list_reserve.csv", row.names = FALSE)
## sp during lockdown
write.csv(list('Species' = lockdown), "01_Analyses_teleo/03_Outputs/Species_list_lockdown.csv", row.names = FALSE)