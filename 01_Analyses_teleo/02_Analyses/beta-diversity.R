library(tidyverse)
library(betapart)
library(ggplot2)
library(reshape2)
library(ggpubr)

####################################################################################
# Beta-diversity between ports

# load data
biodiv_port=read.csv("01_Analyses_teleo/00_data/matrice_teleo_port.csv", row.names=1) %>%
  t(.) %>%
  as.data.frame(.) %>%
  rownames_to_column(var='code_spygen')

species_to_keep_ports <- readRDS("01_Analyses_teleo/00_data/species_to_keep_ports.RDS")

biodiv_port <- biodiv_port[, c("code_spygen", species_to_keep_ports)]

meta_port=read.csv("00_Metadata/metadata_port.csv", sep=";") 


biodiv_port <- left_join(biodiv_port, meta_port[,c("code_spygen", "site")])


# Biodiv per port
biodiv_port <- biodiv_port[,-1]

Site <- unique(biodiv_port$site)
df_port <- data.frame(matrix(0, ncol = 122, nrow = 7))
colnames(df_port) <- colnames(biodiv_port[,-ncol(biodiv_port)])

for (i in 1:length(Site)) {
  df <- biodiv_port %>%
    filter(site == Site[i])
  df_port[i,] <- as.data.frame(t(colSums(df[,-ncol(df)])))
}

rownames(df_port) <- Site
df_port[df_port > 1] <- 1


# Compute beta diversity
b <- betapart.core(df_port)
beta <- beta.pair(b, "jaccard")
Beta_port <- data.frame(matrix(0, ncol=3, nrow=21))
colnames(Beta_port) <- c("Total", "Turnover", "Nestedness")
Beta_port$Total <- as.numeric(beta$beta.jac)
Beta_port$Turnover <- as.numeric(beta$beta.jtu)
Beta_port$Nestedness <- as.numeric(beta$beta.jne)
Beta_port <- melt(Beta_port)
beta_multi_port <- beta.multi(b, "jaccard")

plot_port <- ggplot(Beta_port)+
  geom_violin(aes(x=variable, y=value, fill=variable), show.legend=F)+
  geom_boxplot(aes(x=variable, y=value), width=0.05)+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  xlab("Beta-diversity component")+
  ggtitle("Beta-diversity between ports")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 10))


#################################################################################
# Beta-diversity outside

# load data
biodiv_out=read.csv("01_Analyses_teleo/00_data/biodiv_milieu_naturel.csv", row.names=1) %>%
  t(.) %>%
  as.data.frame(.) %>%
  rownames_to_column(var='code_spygen')
colnames(biodiv_out) <- gsub(" ", "_", colnames(biodiv_out))

species_to_keep_out <- readRDS("01_Analyses_teleo/00_data/species_to_keep_milieuNat.RDS")

biodiv_out <- biodiv_out[, c("code_spygen", species_to_keep_out)]

meta_out=read.csv("00_Metadata/metadata_milieu_naturel.csv", sep=";") 


biodiv_out <- left_join(biodiv_out, meta_out[,c("ï..SPYGEN_code", "Site")], by=c("code_spygen"="ï..SPYGEN_code"))


# Biodiv per site
biodiv_out <- biodiv_out[,-1]

site <- unique(biodiv_out$Site)
df_out <- data.frame(matrix(0, ncol = 247, nrow = 5))
colnames(df_out) <- colnames(biodiv_out[,-ncol(biodiv_out)])

for (i in 1:length(site)) {
  df <- biodiv_out %>%
    filter(Site == site[i])
  df_out[i,] <- as.data.frame(t(colSums(df[,-ncol(df)])))
}

rownames(df_out) <- site
df_out[df_out > 1] <- 1


# Compute beta diversity
b <- betapart.core(df_out)
beta <- beta.pair(b, "jaccard")
Beta_out <- data.frame(matrix(0, ncol=3, nrow=10))
colnames(Beta_out) <- c("Total", "Turnover", "Nestedness")
Beta_out$Total <- as.numeric(beta$beta.jac)
Beta_out$Turnover <- as.numeric(beta$beta.jtu)
Beta_out$Nestedness <- as.numeric(beta$beta.jne)
Beta_out <- melt(Beta_out)
beta_multi_out <- beta.multi(b, "jaccard")

plot_out <- ggplot(Beta_out)+
  geom_violin(aes(x=variable, y=value, fill=variable), show.legend=F)+
  geom_boxplot(aes(x=variable, y=value), width=0.05)+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  xlab("Beta-diversity component")+
  ggtitle("Beta-diversity between natural sites")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 10))




####################################################################################
# Beta-diversity between all sites

# load data
biodiv=read.csv("01_Analyses_teleo/00_data/matrice_teleo_totale.csv", row.names=1) %>%
  t(.) %>%
  as.data.frame(.) %>%
  rownames_to_column(var='code_spygen')

meta=read.csv("00_Metadata/metadata_tot.csv", sep=";") 


biodiv <- left_join(biodiv, meta[,c("code_spygen", "site")])


# Biodiv per port
biodiv <- biodiv[,-1]
biodiv <- biodiv %>%
  filter(!is.na(site))

Site <- unique(biodiv$site)
df_all <- data.frame(matrix(0, ncol = 184, nrow = 13))
colnames(df_all) <- colnames(biodiv[,-ncol(biodiv)])

for (i in 1:length(Site)) {
  df <- biodiv %>%
    filter(site == Site[i])
  df_all[i,] <- as.data.frame(t(colSums(df[,-ncol(df)])))
}

rownames(df_all) <- Site
df_all[df_all > 1] <- 1


# Compute beta diversity
b <- betapart.core(df_all)
beta <- beta.pair(b, "jaccard")
Beta_all <- data.frame(matrix(0, ncol=3, nrow=78))
colnames(Beta_all) <- c("Total", "Turnover", "Nestedness")
Beta_all$Total <- as.numeric(beta$beta.jac)
Beta_all$Turnover <- as.numeric(beta$beta.jtu)
Beta_all$Nestedness <- as.numeric(beta$beta.jne)
Beta_all <- melt(Beta_all)
beta_multi_all <- beta.multi(b, "jaccard")

plot_all <- ggplot(Beta_all)+
  geom_violin(aes(x=variable, y=value, fill=variable), show.legend=F)+
  geom_boxplot(aes(x=variable, y=value), width=0.05)+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  xlab("Beta-diversity component")+
  ggtitle("Beta-diversity between all sites")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 10))


#####################################################################
# Assemble figure

plot_beta <- ggarrange(plot_port, plot_out, plot_all, ncol=3)

ggsave(plot_beta, file="01_Analyses_teleo/03_Outputs/beta-diversity.png", width = 12, height = 6)





# test for differences in turnover
turnover_port <- Beta_port %>%
  filter(variable=="Turnover")
turnover_out <- Beta_out %>%
  filter(variable=="Turnover")
shapiro.test(turnover_port$value)
shapiro.test(turnover_out$value)

wilcox.test(turnover_port$value, turnover_out$value)
t.test(turnover_port$value, turnover_out$value, var.equal=TRUE) 


nestedness_port <- Beta_port %>%
  filter(variable=="Nestedness")
nestedness_out <- Beta_out %>%
  filter(variable=="Nestedness")
shapiro.test(nestedness_port$value)
shapiro.test(nestedness_out$value)

wilcox.test(nestedness_port$value, nestedness_out$value)
t.test(nestedness_port$value, nestedness_out$value, var.equal=TRUE) 
