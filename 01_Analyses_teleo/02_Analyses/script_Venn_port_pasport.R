library("VennDiagram")
#ici nous avons retiré les sp non identifiées
#nous avons également retiré le saumon d'atlantique qui semble très peu probablement présent (contamination ?)
draw.pairwise.venn(159, 95, 86, c("Hors port", "Port"),fill=c("lightblue","yellow"),scaled=T,fontfamily=c("Courier","Courier","Courier"),cat.fontfamily=c("Courier","Courier"), cat.dist=c(-0.05,-0.009), 
                   cat.cex=c(1.4,1.4),cex=c(1.4,1.4,1.4),ext.pos=c(2,10))



