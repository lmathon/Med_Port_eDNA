### Load the libraries
library(dplyr)
library(rsq)
library(fastDummies)
library(tidyr)
library(tidyverse)
library(purrr)
library(scales)
library(margins)
library(ggplot2)
library(viridis)
library(spdep)
library(ggpubr)
library(gridExtra)
library(spaMM)
library(DHARMa)
library(visreg)
library(multcompView)
library(multcomp)
library(MuMIn)
library(visreg)
library(sjPlot)
library(car)
library(MASS)
library(nlme)
library(progressr)
library(gplots)
library(piecewiseSEM)
library(emmeans)


# load data ports
ind_ports <- read.csv("01_Analyses_teleo/00_data/indicators_ports_2022_per_filter.csv", header=T, row.names=1)

# Charger datas milieu naturel (= réserve et hors réserve)
meta_nat <- read.csv2("00_Metadata/metadata_milieu_naturel.csv", header=T)
ind_nat <- read.csv("01_Analyses_teleo/00_data/indicators_milieu_naturel.csv", header=T, row.names = 1) %>%
  rownames_to_column(var="Sample") %>%
  inner_join(meta_nat, by=c("Sample"="SPYGEN_code")) %>%
  rename(habitat = protection,
         X = longitude_start_DD,
         Y = latitude_start_DD) %>%
  column_to_rownames(var="Sample") %>%
  dplyr::select(colnames(ind_ports), "habitat", "Confinement", "X", "Y") %>%
  mutate(Confinement = case_when(Confinement == 'Y' ~ 'Lockdown',
                         TRUE ~ 'Unlock'))

### dataframe ports
meta_ports <- read.csv("00_Metadata/metadata_port.csv", header=T) %>%
  filter(habitat == "Port")
ind_ports <- ind_ports %>%
  rownames_to_column(var="Sample") %>%
  inner_join(meta_ports, by=c("Sample"="code_spygen")) %>%
  mutate(Confinement = "Unlock") %>%
  column_to_rownames(var="Sample") %>%
  dplyr::select(colnames(ind_nat))

# Combiner les deux datasets
ind_all <- rbind(ind_nat,ind_ports) %>%
  mutate(DeBRa = log10(DeBRa)) %>%
  mutate_at('Confinement', as.factor)  %>%
  mutate_at('X', as.numeric)  %>%
  mutate_at('Y', as.numeric) %>%
  mutate(Reserve = case_when(habitat == 'reserve' ~ 'Y',
                                 TRUE ~ 'N')) %>%
  mutate(Port = case_when(habitat == 'Port' ~ 'Y',
                             TRUE ~ 'N')) %>%
  mutate_at('Reserve', as.factor) %>%
  mutate_at('Port', as.factor) 


ind_all$habitat <- factor(ind_all$habitat, levels=c("reserve", "outside", "Port"))
ind_all$Confinement <- factor(ind_all$Confinement, labels=c("Unlock", "Lockdown"))

ind_names <- c("Total species richness", "Cryptobenthic species richness", "Threatened species richness", "Commercial species richness")  


### Matrix of spatial coordinates
coords<- ind_all %>%
  dplyr::select(X, Y) %>%
  as.matrix(.)



###################################################################
## Fit models
###################################################################
# Fit the model
mod <- list()
mod[[1]] <- fitme(R ~ Reserve * Confinement + Port + Matern(1|X+Y),
                  family = "gaussian",  data = ind_all)

mod[[2]] <- fitme(Crypto ~ Reserve * Confinement + Port  + Matern(1|X+Y),
                  family = "gaussian",  data = ind_all)

mod[[3]] <- fitme(RedList ~ Reserve * Confinement + Port  + Matern(1|X+Y),
                  family = "poisson",  data = ind_all)

mod[[4]] <- fitme(Commercial ~ Reserve * Confinement + Port  + Matern(1|X+Y),
                  family = "gaussian",  data = ind_all)


names(mod) <- ind_names
saveRDS(mod, file="01_Analyses_teleo/02_Analyses/models_spaMM_indicators.RDS")




#################################################################################
## Get coefficients and conditional effects of covariates
#################################################################################
mod <- readRDS("01_Analyses_teleo/02_Analyses/models_spaMM_indicators.RDS")


## post-fit evaluation
## create the result dataframe
coef <- matrix(0,5,8, dimnames=list(names(test_all$coef), paste(rep(names(mod),each=2), c("_coef", "_pval"))))
R2 <- vector()


# initiate indices
c=1
p=2

for (i in 1:4) {
  ### test of H_0: all regression coefficients are zero
  test_all <- summary(glht(mod[[i]], coef.=fixef.HLfit))
  coef[,c] <- test_all$coef
  coef[,p] <- test_all$test$pvalues
  
  R2[i] <- cor(ind_all$R,predict(mod[[i]]))^2
  
  c=c+2
  p=p+2
}

names(R2) <- ind_names

write.csv(coef, file="01_Analyses_teleo/02_Analyses/spaMM_indicators_variables_coefficients.csv")
saveRDS(R2, file="01_Analyses_teleo/02_Analyses/spaMM_indicators_R2.RDS")

#############################################################################################
### Partial effects
##########################################################################################
coef <- read.csv("01_Analyses_teleo/02_Analyses/spaMM_indicators_variables_coefficients.csv", header=T, row.names=1)
r2 <- readRDS("01_Analyses_teleo/02_Analyses/spaMM_indicators_R2.RDS")

pdep <- list()

# Extract conditional effects of variables
for (i in 1:4) {
 pdep_port <- pdep_effects(mod[[i]],focal_var="Port")
 pdep_port$focal_var <- c("out", "Port")
 pdep_port$Confinement <- "Unlock"
  
  tmp <- ind_all %>%
    mutate(Confinement = "Unlock")
  
  pdep_1 <- pdep_effects(mod[[i]],focal_var="Reserve", newdata=tmp)
  pdep_1$focal_var <- c("outside", "reserve")
  
  tmp <- ind_all %>%
    mutate(Confinement = "Lockdown")
  
  pdep_2 <- pdep_effects(mod[[i]],focal_var="Reserve", newdata=tmp)
  pdep_2$focal_var <- c("outside", "reserve")
  
  pdep_1$Confinement <- "Lockdown" ; pdep_2$Confinement <- "Unlock"  
  pdep[[i]] <- rbind(pdep_1,pdep_2, pdep_port[2,])
  
  
}

names(pdep) <-  names(mod)
saveRDS(pdep, "01_Analyses_teleo/02_Analyses/spaMM_Conditional_effect_variables.RDS")

#####################################################################################
### Figures 3a 
#####################################################################################
pdep <- readRDS("01_Analyses_teleo/02_Analyses/spaMM_Conditional_effect_variables.RDS")

mytheme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        axis.text.y=element_text(colour="black",size=14),
        axis.text.x=element_text(colour="black",size=15),
        axis.title.y=element_text(colour="black",size=14),
        strip.text = element_text(size = 15),
        axis.title.x=element_blank(),
        plot.title = element_text(colour="black", size=15, face="bold", hjust=0.5),
        plot.subtitle = element_text(colour="black", size=15,  hjust=0.5)) 

# initiate plot list
p <- list()

## Loop to draw plots
for (i in 1:4) {
  v_i <- pdep[[i]]
  v_i$focal_var <- factor(v_i$focal_var, levels=c("reserve", "outside", "Port"))
  v_i$Confinement <- factor(v_i$Confinement, levels=c("Unlock", "Lockdown"))
  R2round <- round(R2[i], digit=3)
  
  p[[i]] <- ggplot(v_i,aes(x=focal_var, y = pointp))+
    #geom_jitter(data = v_i$res, mapping = aes(x = Protection, y = visregRes), colour="grey50", size=0.1) +
    geom_crossbar(aes(ymin=low, ymax=up, fill = focal_var), 
                  linetype=0, alpha=0.3) +
    facet_wrap(~Confinement, scales="free_x") + # argument scales = "free_x" removes the empty factors (no Port during lockdown)
    geom_errorbar(aes(ymin=pointp, ymax=pointp, colour = focal_var), 
                  linetype=1, width=0.9, lwd=1.2) +
    scale_colour_manual(values = c("darkblue", "cyan4",  "darkorange")) +
    scale_fill_manual(values =  c("darkblue", "cyan4",  "darkorange")) +
    labs(y="Species richness",
         title=ind_names[i],
         subtitle = bquote(R^2 == .(R2round))) +
    scale_x_discrete(labels = c("Reserve", "Fished", "Port")) +
    mytheme 
}

## Save plot
png("01_Analyses_teleo/03_Outputs/Figure3a.png", 
    width = 1200, height = 1000, ) 
do.call(grid.arrange,c(p, list(ncol=2)))
dev.off()
