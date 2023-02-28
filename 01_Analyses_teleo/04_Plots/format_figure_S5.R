library(ggpubr)

load("01_Analyses_teleo/04_Plots/WebFig5_beta_diversity.RData")

Fig_S5 <- ggarrange(plot_port, plot_out, plot_all, ncol=3, labels=c("(a)","b)","(c)"))

ggsave(Fig_S5, file="01_Analyses_teleo/04_Plots/Fig_S5.jpeg", width=12, height = 6)
