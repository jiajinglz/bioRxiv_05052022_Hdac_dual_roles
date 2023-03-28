################################################################################
########Plot proportional Hdac1 clusters with AC VG pan-H3Kac enrichment########
################################################################################
dataTable <- data.frame(matrix(nrow = 5, ncol = 4))
dataTable[, 1] <- c("All Clusters", "Cluster I", "Cluster II", "Cluster III", "Cluster IV")
dataTable[, 2] <- c(8100, 1252, 282, 5560, 1006)
dataTable[, 3] <- c(6301, 1611, 50, 4471, 169)
dataTable[, 4] <- c(9041, 685, 1057, 3638, 3661)
colnames(dataTable) <- c("clusters", "LC", "NL", "NP")
#############################################
############make plotable table##############
#############################################
plot_table <- data.frame(matrix(nrow=15, ncol=3))
colnames(plot_table) <- c("type", "value", "cat")
plot_table[1:3, 1] <- c("All Clusters")
plot_table[1:3, 2] <- c(8100, 6301, 9041) ###fill from datatable##
plot_table[1:3, 3] <- c("LC", "NL", "NP")
plot_table[4:6, 1] <- c("Cluster I")
plot_table[4:6, 2] <- c(1251, 1611, 685) ###fill from datatable##
plot_table[4:6, 3] <- c("LC", "NL", "NP")
plot_table[7:9, 1] <- c("Cluster II")
plot_table[7:9, 2] <- c(282, 50, 1057) ###fill from datatable##
plot_table[7:9, 3] <- c("LC", "NL", "NP")
plot_table[10:12, 1] <- c("Cluster III")
plot_table[10:12, 2] <- c(5560, 4471, 3638) ###fill from datatable##
plot_table[10:12, 3] <- c("LC", "NL", "NP")
plot_table[13:15, 1] <- c("Cluster IV")
plot_table[13:15, 2] <- c(1006, 169, 3661) ###fill from datatable##
plot_table[13:15, 3] <- c("LC", "NL", "NP")
##########################################################
plot_table$cat <- factor(plot_table$cat, levels = c("LC", "NL", "NP"))

library(ggplot2)
p <- ggplot(plot_table, aes(x= type, y= value, fill= cat, width= 0.65))+
  geom_bar(position= "fill", stat= "identity")+
  scale_fill_manual(values = c("#87CEEB", "#EEE8AA", "#D3D3D3"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 0.5)+
  theme(axis.text = element_text(color="black"), text = element_text(size=9))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))
p

ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/panH3Kac/differential peak overlap hdac1 clusters/bar plot proportion of AC VG specific pan-H3Kac on hdac1 clusters.tiff', p, device = "tiff", dpi = 200)

















