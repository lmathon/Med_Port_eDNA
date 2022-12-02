### Load the libraries
library(dplyr)
library(rsq)
library(fastDummies)
library(tidyr)
library(purrr)
library(scales)
library(margins)
library(ggplot2)
library(gridExtra)
library(viridis)
library(spdep)

## Load data
data <- read.csv("01_Analyses_teleo/00_data/teleo_presence.csv", header=T, row.names=1) %>%
  filter(rowSums(.) > 0) %>%
  t(.)  %>%
  as.data.frame(.)

meta <- read.csv("00_Metadata/metadata_port.csv", header=T)  %>%
  filter(habitat != "BIOHUT_port")

ind_ports <- read.csv("01_Analyses_teleo/00_data/indicators_ports_2022_per_filter.csv", header=T, row.names=1) %>%
  rownames_to_column(var="code_spygen")

mydata <- ind_ports %>%
  inner_join(meta, by="code_spygen")

### Create a data set of response variables (indicators)
Y <- mydata %>%
  dplyr::select(R,FD,RedList, Commercial) %>%
  # express vqriables as proportion of species richness
  mutate(across(c(FD, RedList,  Commercial), ~ .x/R)) 

### Create a data set of explanatory variables (environment + protection)
X <- mydata %>%
  dplyr::select(port_propre, Campaign,surface_couverte_ha, lineaire_exterieur_m ) %>%
  # transform all character columns to factors
  mutate_if(sapply(., is.character), as.factor)

### Create a vector of site names
sites <- mydata %>%
  dplyr::pull(site)


### Matrix of spatial coordinates
coords<- mydata %>%
  select(X, Y) %>%
  as.matrix(.)

### Create a vector of indicator names
ind_names <- c("Total species richness", "Functional diversity", 
               "Threatened species richness", "Commercial species richness")  


# Create a vector of distributions : Poisson for RedList and Chondri, 
# and Gaussian for all the others
distri <- c("gaussian", "gaussian", 
            "quasipoisson", "gaussian")


#####################################################################################
## Fit the GLMs
#####################################################################################
# create a list to store the models
glm_results <- list()

# Create a vector to store the R-squared
r2 <- rep(NA, ncol(Y), names=ind_names)

# Create a matrix to store the coefficient and pvalues
coef <- matrix(NA, ncol(X), (2*ncol(Y)),
               dimnames= list(colnames(X),
                              c(paste(ind_names,"coef", sep="_" ), paste(ind_names,"pval", sep="_" ))))

# Make a loop to do the glm for each indicator
for (i in 1:ncol(Y)) {
  data_i <- cbind.data.frame(Y[,i], X) # create a dataframe with the data for indicator i
  glm_results[[i]] <- m_i  <- glm(data_i[,1] ~ port_propre + Campaign + surface_couverte_ha + 
                                    lineaire_exterieur_m , 
             data=data_i , family = distri[i], na.action = na.omit )
  
  # Calculate r2 and store it in the result vector
  r2[i] <- rsq(m_i, adj=T)
  
  # Store coefficient in the result matrix
  coef[,i] <- summary(m_i)$coefficients[2:5, 1]
  # Store ANOVA pval in the results
  coef[,(i+ncol(Y))] <- car::Anova(m_i,test.statistic="F")$"Pr(>F)"[1:4]
}


####################################################################
#check Spatial autocorrelation 
lstw <- nb2listw((knn2nb(knearneigh(coords, k=1, longlat = T, 
                                    use_kd_tree=F))))


moran_p <- vector()
for (i in 1:4) {
  m_i <- glm_results[[i]]
  moran_i <-moran.test(m_i$residuals, lstw, 
                       alternative="two.sided")
  moran_p[i] <- moran_i$p.value
}


moran_res <- cbind.data.frame(ind_names,moran_p)
moran_res[which(moran_res$moran_p < 0.05),] # Commercial
mean(moran_p) 
