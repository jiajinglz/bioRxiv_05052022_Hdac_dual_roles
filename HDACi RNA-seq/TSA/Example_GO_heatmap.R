#############################################################
######################Heatmap for GO#########################
#############################################################
list <- read.csv("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 rna seq/ectoderm/AC_GO_AllLists.csv")
table <- data.frame(matrix(nrow=15, ncol=3))
table[, 1] <- c("AC")
go_term <- data.frame(paste(list$GO, list$Description))
table[, 2] <- go_term[1:15, ]
table[, 3] <- list[1:15, 6]
colnames(table) <- c("region", "term", "LogP")
library(ggplot2)
table$term <- factor(table$term, levels=rev(table$term))
p <- ggplot(table, aes(region, term, fill= LogP)) + 
     geom_tile(color = "black", lwd = 0.5, linetype = 1)+
     scale_fill_distiller(palette = "RdPu")+
     coord_fixed()+
     scale_y_discrete(position = "right")+
     theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.ticks.y = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           panel.background = element_blank(),
           axis.text.y = element_text(size = 30, color = "black", margin = margin(-20, -20, -20, -20)))
p

ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 rna seq/ectoderm/AC_GO_heatmap.tiff', p, device = "tiff", dpi = 200, width = 40, height = 20, units = c("cm"))








