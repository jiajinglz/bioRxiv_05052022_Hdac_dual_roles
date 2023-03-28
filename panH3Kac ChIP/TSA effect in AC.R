######################################################################
######################################################################
####################TSA affect panH3Kac in AC#########################
######################################################################
ac_dmso_all <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_HDAC1_combined_peak_1kb.AC_DMSOcounts", header=FALSE)
ac_tsa_all <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_HDAC1_combined_peak_1kb.AC_TSAcounts", header=FALSE)

ac_dmso_c1 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_hdac1_h3k27ac_h3k27me3_overlap_peak_1kb.AC_DMSOcounts", header=FALSE)
ac_tsa_c1 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_hdac1_h3k27ac_h3k27me3_overlap_peak_1kb.AC_TSAcounts", header=FALSE)

ac_dmso_c2 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_hdac1_only_h3k27me3_peak_1kb.AC_DMSOcounts", header=FALSE)
ac_tsa_c2 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_hdac1_only_h3k27me3_peak_1kb.AC_TSAcounts", header=FALSE)

ac_dmso_c3 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_hdac1_only_h3k27ac_peak_1kb.AC_DMSOcounts", header=FALSE)
ac_tsa_c3 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_hdac1_only_h3k27ac_peak_1kb.AC_TSAcounts", header=FALSE)

ac_dmso_c4 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_hdac1_neither_h3k27ac_h3k27me3_peak_1kb.AC_DMSOcounts", header=FALSE)
ac_tsa_c4 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_hdac1_neither_h3k27ac_h3k27me3_peak_1kb.AC_TSAcounts", header=FALSE)

ac_dmso_c5 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_random_regions_1kb.AC_DMSOcounts", header=FALSE)
ac_tsa_c5 <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/st9_st105_random_regions_1kb.AC_TSAcounts", header=FALSE)
############################################################################################
##########################effect by TSA in each cluster#####################################
############################################################################################
rownames(ac_dmso_all) <- paste(ac_dmso_all$V1, ac_dmso_all$V2, ac_dmso_all$V3, sep = "\t")
ac_dmso_all <- ac_dmso_all[, 4, drop = F]
ac_dmso_all[, 1] <- ac_dmso_all$V4/3.86

rownames(ac_dmso_c1) <- paste(ac_dmso_c1$V1, ac_dmso_c1$V2, ac_dmso_c1$V3, sep = "\t")
ac_dmso_c1 <- ac_dmso_c1[, 4, drop = F]
ac_dmso_c1[, 1] <- ac_dmso_c1$V4/3.86
mean(ac_dmso_c1$V4)

rownames(ac_dmso_c2) <- paste(ac_dmso_c2$V1, ac_dmso_c2$V2, ac_dmso_c2$V3, sep = "\t")
ac_dmso_c2 <- ac_dmso_c2[, 4, drop = F]
ac_dmso_c2[, 1] <- ac_dmso_c2$V4/3.86
mean(ac_dmso_c2$V4)

rownames(ac_dmso_c3) <- paste(ac_dmso_c3$V1, ac_dmso_c3$V2, ac_dmso_c3$V3, sep = "\t")
ac_dmso_c3 <- ac_dmso_c3[, 4, drop = F]
ac_dmso_c3[, 1] <- ac_dmso_c3$V4/3.86
mean(ac_dmso_c3$V4)

rownames(ac_dmso_c4) <- paste(ac_dmso_c4$V1, ac_dmso_c4$V2, ac_dmso_c4$V3, sep = "\t")
ac_dmso_c4 <- ac_dmso_c4[, 4, drop = F]
ac_dmso_c4[, 1] <- ac_dmso_c4$V4/3.86
mean(ac_dmso_c4$V4)

rownames(ac_dmso_c5) <- paste(ac_dmso_c5$V1, ac_dmso_c5$V2, ac_dmso_c5$V3, sep = "\t")
ac_dmso_c5 <- ac_dmso_c5[, 7, drop = F]
ac_dmso_c5[, 1] <- ac_dmso_c5$V7/3.86

library(dplyr)
ac_dmso_all <- ac_dmso_all %>% mutate(type=c("all"))
ac_dmso_c1 <- ac_dmso_c1 %>% mutate(type=c("Cluster1"))
ac_dmso_c2 <- ac_dmso_c2 %>% mutate(type=c("Cluster2"))
ac_dmso_c3 <- ac_dmso_c3 %>% mutate(type=c("Cluster3"))
ac_dmso_c4 <- ac_dmso_c4 %>% mutate(type=c("Cluster4"))
ac_dmso_c5 <- ac_dmso_c5 %>% mutate(type=c("Ran_regions"))
colnames(ac_dmso_c5) <- c("V4", "type")

ac_dmso_vpTable <- rbind(ac_dmso_all, ac_dmso_c1, ac_dmso_c2, 
                         ac_dmso_c3, ac_dmso_c4, ac_dmso_c5)
colnames(ac_dmso_vpTable) <- c("value", "type")

##########################

rownames(ac_tsa_all) <- paste(ac_tsa_all$V1, ac_tsa_all$V2, ac_tsa_all$V3, sep = "\t")
ac_tsa_all <- ac_tsa_all[, 4, drop = F]
ac_tsa_all[, 1] <- ac_tsa_all$V4/1.34

rownames(ac_tsa_c1) <- paste(ac_tsa_c1$V1, ac_tsa_c1$V2, ac_tsa_c1$V3, sep = "\t")
ac_tsa_c1 <- ac_tsa_c1[, 4, drop = F]
ac_tsa_c1[, 1] <- ac_tsa_c1$V4/1.34

rownames(ac_tsa_c2) <- paste(ac_tsa_c2$V1, ac_tsa_c2$V2, ac_tsa_c2$V3, sep = "\t")
ac_tsa_c2 <- ac_tsa_c2[, 4, drop = F]
ac_tsa_c2[, 1] <- ac_tsa_c2$V4/1.34

rownames(ac_tsa_c3) <- paste(ac_tsa_c3$V1, ac_tsa_c3$V2, ac_tsa_c3$V3, sep = "\t")
ac_tsa_c3 <- ac_tsa_c3[, 4, drop = F]
ac_tsa_c3[, 1] <- ac_tsa_c3$V4/1.34

rownames(ac_tsa_c4) <- paste(ac_tsa_c4$V1, ac_tsa_c4$V2, ac_tsa_c4$V3, sep = "\t")
ac_tsa_c4 <- ac_tsa_c4[, 4, drop = F]
ac_tsa_c4[, 1] <- ac_tsa_c4$V4/1.34

rownames(ac_tsa_c5) <- paste(ac_tsa_c5$V1, ac_tsa_c5$V2, ac_tsa_c5$V3, sep = "\t")
ac_tsa_c5 <- ac_tsa_c5[, 7, drop = F]
ac_tsa_c5[, 1] <- ac_tsa_c5$V7/1.34


library(dplyr)
ac_tsa_all <- ac_tsa_all %>% mutate(type=c("all"))
ac_tsa_c1 <- ac_tsa_c1 %>% mutate(type=c("Cluster1"))
ac_tsa_c2 <- ac_tsa_c2 %>% mutate(type=c("Cluster2"))
ac_tsa_c3 <- ac_tsa_c3 %>% mutate(type=c("Cluster3"))
ac_tsa_c4 <- ac_tsa_c4 %>% mutate(type=c("Cluster4"))
ac_tsa_c5 <- ac_tsa_c5 %>% mutate(type=c("Ran_regions"))
colnames(ac_tsa_c5) <- c("V4", "type")

ac_tsa_vpTable <- rbind(ac_tsa_all, ac_tsa_c1, ac_tsa_c2, 
                        ac_tsa_c3, ac_tsa_c4, ac_tsa_c5)
colnames(ac_tsa_vpTable) <- c("value", "type")

###############################################################
############Plot normalized reads in DMSO and TSA##############
###############################################################
ac_dmso_vpTable <- ac_dmso_vpTable %>% mutate(condition=c("DMSO"))
ac_tsa_vpTable <- ac_tsa_vpTable %>% mutate(condition=c("TSA"))
ac_vpTable <- rbind(ac_dmso_vpTable, ac_tsa_vpTable)

library(ggplot2)
library(RColorBrewer)
p <- ggplot(ac_vpTable, aes(x= type, y= log(value+1, 2)))+ 
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

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/pan-H3Kac_AC_hdac1_clusters.tiff', p, device = "tiff", dpi = 200, width = 20, height = 10, units = c("cm"))

###############################################################
######Calculate the exact difference between DMSO and TSA######
###############################################################
dif_all <- data.frame(matrix(nrow=nrow(ac_dmso_all), ncol=2))
rownames(dif_all) <- rownames(ac_dmso_all)
dif_all[, 1] <- ac_dmso_all$type
dif_all[, 2] <- ac_tsa_all$V4 - ac_dmso_all$V4
colnames(dif_all) <- c("type", "dif_value")
min(dif_all$dif_value)

dif_c1 <- data.frame(matrix(nrow=nrow(ac_dmso_c1), ncol=2))
rownames(dif_c1) <- rownames(ac_dmso_c1)
dif_c1[, 1] <- ac_dmso_c1$type
dif_c1[, 2] <- ac_tsa_c1$V4 - ac_dmso_c1$V4
colnames(dif_c1) <- c("type", "dif_value")

dif_c2 <- data.frame(matrix(nrow=nrow(ac_dmso_c2), ncol=2))
rownames(dif_c2) <- rownames(ac_dmso_c2)
dif_c2[, 1] <- ac_dmso_c2$type
dif_c2[, 2] <- ac_tsa_c2$V4 - ac_dmso_c2$V4
colnames(dif_c2) <- c("type", "dif_value")

dif_c3 <- data.frame(matrix(nrow=nrow(ac_dmso_c3), ncol=2))
rownames(dif_c3) <- rownames(ac_dmso_c3)
dif_c3[, 1] <- ac_dmso_c3$type
dif_c3[, 2] <- ac_tsa_c3$V4 - ac_dmso_c3$V4
colnames(dif_c3) <- c("type", "dif_value")

dif_c4 <- data.frame(matrix(nrow=nrow(ac_dmso_c4), ncol=2))
rownames(dif_c4) <- rownames(ac_dmso_c4)
dif_c4[, 1] <- ac_dmso_c4$type
dif_c4[, 2] <- ac_tsa_c4$V4 - ac_dmso_c4$V4
colnames(dif_c4) <- c("type", "dif_value")

dif_c5 <- data.frame(matrix(nrow=nrow(ac_dmso_c5), ncol=2))
rownames(dif_c5) <- rownames(ac_dmso_c5)
dif_c5[, 1] <- ac_dmso_c5$type
dif_c5[, 2] <- ac_tsa_c5$V4 - ac_dmso_c5$V4
colnames(dif_c5) <- c("type", "dif_value")

dif_table <- rbind(dif_all, dif_c1, dif_c2, dif_c3, dif_c4, dif_c5)

library(ggplot2)
library(RColorBrewer)
p <- ggplot(dif_table, aes(x= type, y= log(dif_value + 1,2), fill= type))+ 
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

mean(dif_c1$dif_value)
summary(dif_c1$dif_value)
mean(dif_c2$dif_value)
summary(dif_c2$dif_value)
mean(dif_c3$dif_value)
summary(dif_c3$dif_value)
mean(dif_c4$dif_value)
summary(dif_c4$dif_value)
mean(dif_c5$dif_value)
summary(dif_c5$dif_value)

t.test(dif_c2$dif_value, dif_c3$dif_value)

summary(subset(dif_c5, dif_c5$dif_value< 1))

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/pan-H3Kac_VG_hdac1_clusters_direct_subtraction.tiff', p, device = "tiff", dpi = 200, width = 20, height = 10, units = c("cm"))

#############################################################
####################foldChange of TSA/DMSO###################
#############################################################
foldChange <- data.frame(matrix(nrow= nrow(ac_dmso_vpTable), ncol = 2))
rownames(foldChange) <- row.names(ac_dmso_vpTable)
foldChange[, 1] <- (ac_tsa_vpTable$value/ (ac_dmso_vpTable$value+1))
foldChange[, 2] <- (ac_tsa_vpTable$type)
colnames(foldChange) <- c("fold", "type")

summary(subset(foldChange, foldChange$type == "Cluster1"))
summary(subset(foldChange, foldChange$type == "Cluster2"))
summary(subset(foldChange, foldChange$type == "Cluster3"))
summary(subset(foldChange, foldChange$type == "Cluster4"))
summary(subset(foldChange, foldChange$type == "Ran_regions"))
summary(subset(foldChange, foldChange$type == "all"))


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

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/ac_counts/pan-H3Kac_AC_hdac1_clusters_FoldIncrease_by_TSA_Scaled.tiff', p, device = "tiff", dpi = 200, width = 20, height = 10, units = c("cm"))

x<- subset(foldChange, foldChange$type== "Ran_regions")
x<- x[complete.cases(x),]
x<- x[is.finite(x$fold),]
y<- subset(foldChange, foldChange$type== "all")
y<- y[complete.cases(y),]
y<- y[is.finite(y$fold),]
t.test(x$fold, y$fold, alternative = "two.sided")



