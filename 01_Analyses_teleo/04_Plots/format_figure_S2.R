library(ggpubr)

load("01_Analyses_teleo/04_Plots/WebFig2_a_families_count.RData")
load("01_Analyses_teleo/04_Plots/WebFig2_b_families_count_reduced.RData")

Fig_S2 <- ggarrange(plot_family_count, plot_family_count_reduced, 
                    ncol=2, labels=c("(a)","(b)"), common.legend = T, legend = "right")


ggsave(Fig_S2, file="01_Analyses_teleo/04_Plots/Fig_S2.jpeg", width=7, height = 14)
