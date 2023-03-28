######################################################################
######################################################################
####################TSA affect panH3Kac in AC#########################
######################################################################
vg_dmso_all <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_HDAC1_combined_peak_1kb.VG_DMSOcounts", header=FALSE)
vg_tsa_all <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_HDAC1_combined_peak_1kb.VG_TSAcounts", header=FALSE)

vg_dmso_c1 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_hdac1_h3k27ac_h3k27me3_overlap_peak_1kb.VG_DMSOcounts", header=FALSE)
vg_tsa_c1 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_hdac1_h3k27ac_h3k27me3_overlap_peak_1kb.VG_TSAcounts", header=FALSE)

vg_dmso_c2 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_hdac1_only_h3k27me3_peak_1kb.VG_DMSOcounts", header=FALSE)
vg_tsa_c2 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_hdac1_only_h3k27me3_peak_1kb.VG_TSAcounts", header=FALSE)

vg_dmso_c3 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_hdac1_only_h3k27ac_peak_1kb.VG_DMSOcounts", header=FALSE)
vg_tsa_c3 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_hdac1_only_h3k27ac_peak_1kb.VG_TSAcounts", header=FALSE)

vg_dmso_c4 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_hdac1_neither_h3k27ac_h3k27me3_peak_1kb.VG_DMSOcounts", header=FALSE)
vg_tsa_c4 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_hdac1_neither_h3k27ac_h3k27me3_peak_1kb.VG_TSAcounts", header=FALSE)

vg_dmso_c5 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_random_regions_1kb.VG_DMSOcounts", header=FALSE)
vg_tsa_c5 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/st9_st105_random_regions_1kb.VG_TSAcounts", header=FALSE)
############################################################################################
##########################effect by TSA in each cluster#####################################
############################################################################################
rownames(vg_dmso_all) <- paste(vg_dmso_all$V1, vg_dmso_all$V2, vg_dmso_all$V3, sep = "\t")
vg_dmso_all <- vg_dmso_all[, 4, drop = F]
vg_dmso_all[, 1] <- vg_dmso_all$V4/3.95

rownames(vg_dmso_c1) <- paste(vg_dmso_c1$V1, vg_dmso_c1$V2, vg_dmso_c1$V3, sep = "\t")
vg_dmso_c1 <- vg_dmso_c1[, 4, drop = F]
vg_dmso_c1[, 1] <- vg_dmso_c1$V4/3.95
mean(vg_dmso_c1$V4)

rownames(vg_dmso_c2) <- paste(vg_dmso_c2$V1, vg_dmso_c2$V2, vg_dmso_c2$V3, sep = "\t")
vg_dmso_c2 <- vg_dmso_c2[, 4, drop = F]
vg_dmso_c2[, 1] <- vg_dmso_c2$V4/3.95
mean(vg_dmso_c2$V4)

rownames(vg_dmso_c3) <- paste(vg_dmso_c3$V1, vg_dmso_c3$V2, vg_dmso_c3$V3, sep = "\t")
vg_dmso_c3 <- vg_dmso_c3[, 4, drop = F]
vg_dmso_c3[, 1] <- vg_dmso_c3$V4/3.95
mean(vg_dmso_c3$V4)

rownames(vg_dmso_c4) <- paste(vg_dmso_c4$V1, vg_dmso_c4$V2, vg_dmso_c4$V3, sep = "\t")
vg_dmso_c4 <- vg_dmso_c4[, 4, drop = F]
vg_dmso_c4[, 1] <- vg_dmso_c4$V4/3.95
mean(vg_dmso_c4$V4)

rownames(vg_dmso_c5) <- paste(vg_dmso_c5$V1, vg_dmso_c5$V2, vg_dmso_c5$V3, sep = "\t")
vg_dmso_c5 <- vg_dmso_c5[, 7, drop = F]
vg_dmso_c5[, 1] <- vg_dmso_c5$V7/3.95

library(dplyr)
vg_dmso_all <- vg_dmso_all %>% mutate(type=c("all"))
vg_dmso_c1 <- vg_dmso_c1 %>% mutate(type=c("Cluster1"))
vg_dmso_c2 <- vg_dmso_c2 %>% mutate(type=c("Cluster2"))
vg_dmso_c3 <- vg_dmso_c3 %>% mutate(type=c("Cluster3"))
vg_dmso_c4 <- vg_dmso_c4 %>% mutate(type=c("Cluster4"))
vg_dmso_c5 <- vg_dmso_c5 %>% mutate(type=c("Ran_regions"))
colnames(vg_dmso_c5) <- c("V4", "type")

vg_dmso_vpTable <- rbind(vg_dmso_all, vg_dmso_c1, vg_dmso_c2, 
                         vg_dmso_c3, vg_dmso_c4, vg_dmso_c5)
colnames(vg_dmso_vpTable) <- c("value", "type")

##########################

rownames(vg_tsa_all) <- paste(vg_tsa_all$V1, vg_tsa_all$V2, vg_tsa_all$V3, sep = "\t")
vg_tsa_all <- vg_tsa_all[, 4, drop = F]
vg_tsa_all[, 1] <- vg_tsa_all$V4/1.13

rownames(vg_tsa_c1) <- paste(vg_tsa_c1$V1, vg_tsa_c1$V2, vg_tsa_c1$V3, sep = "\t")
vg_tsa_c1 <- vg_tsa_c1[, 4, drop = F]
vg_tsa_c1[, 1] <- vg_tsa_c1$V4/1.13

rownames(vg_tsa_c2) <- paste(vg_tsa_c2$V1, vg_tsa_c2$V2, vg_tsa_c2$V3, sep = "\t")
vg_tsa_c2 <- vg_tsa_c2[, 4, drop = F]
vg_tsa_c2[, 1] <- vg_tsa_c2$V4/1.13

rownames(vg_tsa_c3) <- paste(vg_tsa_c3$V1, vg_tsa_c3$V2, vg_tsa_c3$V3, sep = "\t")
vg_tsa_c3 <- vg_tsa_c3[, 4, drop = F]
vg_tsa_c3[, 1] <- vg_tsa_c3$V4/1.13

rownames(vg_tsa_c4) <- paste(vg_tsa_c4$V1, vg_tsa_c4$V2, vg_tsa_c4$V3, sep = "\t")
vg_tsa_c4 <- vg_tsa_c4[, 4, drop = F]
vg_tsa_c4[, 1] <- vg_tsa_c4$V4/1.13

rownames(vg_tsa_c5) <- paste(vg_tsa_c5$V1, vg_tsa_c5$V2, vg_tsa_c5$V3, sep = "\t")
vg_tsa_c5 <- vg_tsa_c5[, 7, drop = F]
vg_tsa_c5[, 1] <- vg_tsa_c5$V7/1.13


library(dplyr)
vg_tsa_all <- vg_tsa_all %>% mutate(type=c("all"))
vg_tsa_c1 <- vg_tsa_c1 %>% mutate(type=c("Cluster1"))
vg_tsa_c2 <- vg_tsa_c2 %>% mutate(type=c("Cluster2"))
vg_tsa_c3 <- vg_tsa_c3 %>% mutate(type=c("Cluster3"))
vg_tsa_c4 <- vg_tsa_c4 %>% mutate(type=c("Cluster4"))
vg_tsa_c5 <- vg_tsa_c5 %>% mutate(type=c("Ran_regions"))
colnames(vg_tsa_c5) <- c("V4", "type")

vg_tsa_vpTable <- rbind(vg_tsa_all, vg_tsa_c1, vg_tsa_c2, 
                        vg_tsa_c3, vg_tsa_c4, vg_tsa_c5)
colnames(vg_tsa_vpTable) <- c("value", "type")

###############################################################
############Plot normalized reads in DMSO and TSA##############
###############################################################

vg_dmso_vpTable <- vg_dmso_vpTable %>% mutate(condition=c("DMSO"))
vg_tsa_vpTable <- vg_tsa_vpTable %>% mutate(condition=c("TSA"))
vg_vpTable <- rbind(vg_dmso_vpTable, vg_tsa_vpTable)
vg_vpTable [, 1] <- log(vg_vpTable$value+1,2) ## log transformation

library(ggplot2)
library(RColorBrewer)
p <- ggplot(vg_vpTable, aes(x= type, y= value))+ 
  geom_violin(aes(fill= condition), position=position_dodge(.9), trim=T, scale = "area", color= "white")+
  stat_summary(fun=mean, aes(group= condition), position=position_dodge(.9), geom="point", size=2, color= "black")+
  scale_fill_manual(values=c("#4682B4", "#DA70D6"))+
  geom_boxplot(width=0.1, aes(fill= condition), position=position_dodge(.9), color="grey", alpha=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 1/2)
p

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/pan-H3Kac_VG_hdac1_clusters.tiff', p, device = "tiff", dpi = 200, width = 20, height = 10, units = c("cm"))


###############################################################
######Calculate the exact difference between DMSO and TSA######
###############################################################
table <- data.frame(matrix(nrow=nrow(vg_dmso_vpTable), ncol=2))
rownames(table) <- rownames(vg_dmso_vpTable)
table[, 1] <- vg_dmso_vpTable$type
table[, 2] <- vg_tsa_vpTable$value - vg_dmso_vpTable$value
colnames(table) <- c("type", "dif_value")
table[, 2] <- log(table$dif_value, 2)

library(ggplot2)
library(RColorBrewer)
p <- ggplot(table, aes(x= type, y= dif_value, fill= type))+ 
  geom_violin(trim=T, scale = "area", color= "white")+
  stat_summary(fun=mean, geom="point", size=2, color= "black")+
  scale_fill_manual(values=c("#008080", "#008080", "#008080", "#008080", "#008080", "#008080"))+
  geom_boxplot(width=0.1, color="grey", alpha=0.5)+
  geom_hline(yintercept = 0, linetype='dotted', col = 'red')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 1/2)
p

x <- subset(table, table$type == "Cluster4")
mean(x$dif_value)
y <- subset(table, table$type == "Ran_regions")
mean(y$dif_value)

t.test(x$dif_value, y$dif_value, alternative= "two.sided")

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/pan-H3Kac_AC_hdac1_clusters_direct_subtraction.tiff', p, device = "tiff", dpi = 200, width = 20, height = 10, units = c("cm"))


#############################################################
####################foldChange of TSA/DMSO###################
#############################################################

foldChange <- data.frame(matrix(nrow= nrow(vg_dmso_vpTable), ncol = 2))
rownames(foldChange) <- row.names(vg_dmso_vpTable)
foldChange[, 1] <- (vg_tsa_vpTable$value/ (vg_dmso_vpTable$value+1))
foldChange[, 2] <- (vg_tsa_vpTable$type)
colnames(foldChange) <- c("fold", "type")

library(ggplot2)
library(RColorBrewer)
p <- ggplot(foldChange, aes(x= type, y= log(fold, 2), fill= type))+ 
  geom_violin(trim=T, scale = "area", color= "white")+
  stat_summary(fun=mean, geom="point", size=2, color= "black")+
  scale_fill_manual(values=c("burlywood1", "burlywood1", "burlywood1", "burlywood1", "burlywood1", "burlywood1"))+
  geom_boxplot(width=0.1, color="grey", alpha=0.5)+
  geom_hline(yintercept = 0, linetype='dotted', col = 'red')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 1/2)
p

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/vg_counts/pan-H3Kac_VG_ChIP_4_hdac1_clusters_FoldIncrease.tiff', p, device = "tiff", dpi = 200, width = 20, height = 10, units = c("cm"))

x<- subset(foldChange, foldChange$type== "all")
x<- x[complete.cases(x),]
x<- x[is.finite(x$fold),]
y<- subset(foldChange, foldChange$type== "Ran_regions")
y<- y[complete.cases(y),]
y<- y[is.finite(y$fold),]
t.test(x$fold, y$fold, alternative = "two.sided")






