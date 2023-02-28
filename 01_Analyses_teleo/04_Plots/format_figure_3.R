## ---------------------------
##
## Format figure 3
##
## Date Created: 2023-02-28
##
## ---------------------------

#_libs
library(ggpubr)

# Load RData
load('01_Analyses_teleo/04_Plots/Fig3_ab_venn.RData')
load('01_Analyses_teleo/04_Plots/Fig3_cf.RData')
load('01_Analyses_teleo/04_Plots/Fig3_gh_dbRDA_totale.RData')

# Format Figure 3
venn_a <- ggpubr::ggarrange(venn_portreservfish, ncol = 1)
venn_b <- ggpubr::ggarrange(venn_portlockun, ncol = 1)
fig3 <- ggpubr::ggarrange(venn_a, venn_b, p[[1]], p[[2]], p[[3]], p[[4]], grda_sites_total, grda_sites_turnover, 
                          ncol = 2, nrow = 4, labels = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"))

# Save: tiff, jpeg

tiff(filename = "01_Analyses_teleo/04_Plots/Fig3.tiff", height = 5000, width = 3200, res = 300, units = "px")
fig3
dev.off()

jpeg(filename = "01_Analyses_teleo/04_Plots/Fig3.jpeg", height = 5000, width = 3200, res = 300, units = "px", quality = 100)
fig3
dev.off()
