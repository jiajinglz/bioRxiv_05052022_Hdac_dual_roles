################################################################################
###############Hdac1 Cluster 1 TSA effect AC vs VG##############################
################################################################################
c1_AC_DMSO <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/hdac1_cluster1/st9_st105_hdac1_h3k27ac_h3k27me3_overlap_peak_1kb.AC_DMSOcounts", header=FALSE)
c1_AC_TSA <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/hdac1_cluster1/st9_st105_hdac1_h3k27ac_h3k27me3_overlap_peak_1kb.AC_TSAcounts", header=FALSE)
c1_VG_DMSO <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/hdac1_cluster1/st9_st105_hdac1_h3k27ac_h3k27me3_overlap_peak_1kb.VG_DMSOcounts", header=FALSE)
c1_VG_TSA <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/hdac1_cluster1/st9_st105_hdac1_h3k27ac_h3k27me3_overlap_peak_1kb.VG_TSAcounts", header=FALSE)

counts_table <- data.frame(matrix(nrow = nrow(c1_AC_DMSO), ncol = 4))
rownames(counts_table) <- paste(c1_AC_DMSO$V1, c1_AC_DMSO$V2, c1_AC_DMSO$V3, sep = "\t")
counts_table[, 1] <- c1_AC_DMSO$V4/3.86  ##divide by scaling factor from spike-in####
counts_table[, 2] <- c1_AC_TSA$V4/1.34   ##divide by scaling factor from spike-in####
counts_table[, 3] <- c1_VG_DMSO$V4/3.95  ##divide by scaling factor from spike-in####
counts_table[, 4] <- c1_VG_TSA$V4/1.13   ##divide by scaling factor from spike-in####
colnames(counts_table) <- c("AC_DMSO", "AC_TSA", "VG_DMSO", "VG_TSA")


##################################################
###try plot TSA effect based on AC_signal in AC###
##################################################
plot_table <- data.frame(matrix(nrow = nrow(counts_table), ncol = 3))
plot_table[, 1] <- counts_table$AC_DMSO
plot_table[, 2] <- counts_table$AC_TSA/ counts_table$AC_DMSO
plot_table[, 3] <- counts_table$VG_TSA/ counts_table$VG_DMSO
colnames(plot_table) <- c("ac_signal", "AC_fc", "VG_fc")
plot_table <- plot_table[order(plot_table$ac_signal), ]

library(dplyr)
plot_table <- plot_table %>% mutate(1: 3548)
colnames(plot_table) <- c("ac_signal", "AC_fc", "VG_fc", "rank")

library(ggplot2)
library(RColorBrewer)
g <- ggplot(plot_table, aes(x=rank, y=log(AC_fc,2)))+
  geom_point(alpha=0.5, size= 0.5, color= "tomato1")+
  geom_smooth(method = lm, color= "grey50")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(aspect.ratio=1)
g

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/hdac1_cluster1/panH3Kac_AC_foldChange_along_AC_signal.tiff', g, device = "tiff", dpi = 200)



###try plot TSA effect based on VG_signal in VG###
plot_table <- data.frame(matrix(nrow = nrow(counts_table), ncol = 3))
plot_table[, 1] <- counts_table$VG_DMSO
plot_table[, 2] <- counts_table$AC_TSA/ counts_table$AC_DMSO
plot_table[, 3] <- counts_table$VG_TSA/ counts_table$VG_DMSO
colnames(plot_table) <- c("vg_signal", "AC_fc", "VG_fc")
plot_table <- plot_table[order(plot_table$vg_signal), ]

library(dplyr)
plot_table <- plot_table %>% mutate(1: 3548)
colnames(plot_table) <- c("vg_signal", "AC_fc", "VG_fc", "rank")

library(ggplot2)
library(RColorBrewer)
g <- ggplot(plot_table, aes(x=rank, y= log(VG_fc, 2)))+
  geom_point(alpha=0.5, size= 0.5, color= "palegreen3")+
  geom_smooth(method = lm, color= "grey50")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(aspect.ratio=1)
g

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/hdac1_cluster1/panH3Kac_VG_foldChange_along_VG_signal.tiff.tiff', g, device = "tiff", dpi = 200)


######################################################################
###########try to combine previous two plots together#################
############################Violin plot###############################

plot_table <- data.frame(matrix(nrow = nrow(counts_table), ncol = 3))
plot_table[, 1] <- counts_table$VG_DMSO/ counts_table$AC_DMSO ##larger number= VG CRMs
plot_table[, 2] <- counts_table$AC_TSA/ counts_table$AC_DMSO
plot_table[, 3] <- counts_table$VG_TSA/ counts_table$VG_DMSO
colnames(plot_table) <- c("FC", "AC_fc", "VG_fc")
plot_table <- plot_table[order(plot_table$FC), ]
nrow(subset(plot_table, plot_table$FC <= 1/2)) ## number of AC CRMs
nrow(subset(plot_table, plot_table$FC >= 2/1)) ## number of VG CRMs

library(dplyr)
ac_crm <- subset(plot_table, plot_table$FC <= 1/2)
vg_crm <- subset(plot_table, plot_table$FC >= 2/1)
np_crm <- subset(plot_table, plot_table$FC > 1/2 & plot_table$FC < 2/1)

ac_crm <- ac_crm %>% mutate(type= c("AC CRMs"))
vg_crm <- vg_crm %>% mutate(type= c("VG CRMs"))
np_crm <- np_crm %>% mutate(type= c("Ubiquitous CRMs"))

col <- c(2, 4)
ac_crm_acChange <- ac_crm[, col]
ac_crm_acChange <- ac_crm_acChange %>% mutate(condition= c("FC in AC"))
colnames(ac_crm_acChange) <- c("value", "type", "condition")

col <- c(3, 4)
ac_crm_vgChange <- ac_crm[, col]
ac_crm_vgChange <- ac_crm_vgChange %>% mutate(condition= c("FC in VG"))
colnames(ac_crm_vgChange) <- c("value", "type", "condition")

col <- c(2,4)
vg_crm_acChange <- vg_crm[, col]
vg_crm_acChange <- vg_crm_acChange %>% mutate(condition= c("FC in AC"))
colnames(vg_crm_acChange) <- c("value", "type", "condition")

col <- c(3, 4)
vg_crm_vgChange <- vg_crm[, col]
vg_crm_vgChange <- vg_crm_vgChange %>% mutate(condition= c("FC in VG"))
colnames(vg_crm_vgChange) <- c("value", "type", "condition")

col <- c(2, 4)
np_crm_acChange <- np_crm[, col]
np_crm_acChange <- np_crm_acChange %>% mutate(condition= c("FC in AC"))
colnames(np_crm_acChange) <- c("value", "type", "condition")

col <- c(3, 4)
np_crm_vgChange <- np_crm[, col]
np_crm_vgChange <- np_crm_vgChange %>% mutate(condition= c("FC in VG"))
colnames(np_crm_vgChange) <- c("value", "type", "condition")

bp_table <- rbind(ac_crm_acChange, ac_crm_vgChange, 
                  vg_crm_acChange, vg_crm_vgChange,
                  np_crm_acChange, np_crm_vgChange)

library(ggplot2)
library(RColorBrewer)
p <- ggplot(bp_table, aes(x= type, y= log(value,2)))+ 
  geom_violin(aes(fill= condition), position=position_dodge(.9), trim=T, scale = "area", color= "white")+
  stat_summary(fun=mean, aes(group= condition), position=position_dodge(.9), geom="point", size=2, color= "black")+
  scale_fill_manual(values = alpha(c("tomato1", "palegreen3"), 0.7))+
  geom_boxplot(width=0.1, aes(fill= condition), position=position_dodge(.9), color="grey", alpha=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(aspect.ratio = 1/2)

p

#ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/counts by TSA/hdac1_cluster1/vplot_panH3Kac_foldChange_of_spatial_CRMs.tiff', p, device = "tiff", dpi = 200)

a<- subset(bp_table, bp_table$type== "AC CRMs"& bp_table$condition== "FC in AC")
b<- subset(bp_table, bp_table$type== "AC CRMs"& bp_table$condition== "FC in VG")

c<- subset(bp_table, bp_table$type== "Ubiquitous CRMs"& bp_table$condition== "FC in AC")
d<- subset(bp_table, bp_table$type== "Ubiquitous CRMs"& bp_table$condition== "FC in VG")

e<- subset(bp_table, bp_table$type== "VG CRMs"& bp_table$condition== "FC in AC")
f<- subset(bp_table, bp_table$type== "VG CRMs"& bp_table$condition== "FC in VG")

library(lsr)
cohensD(a$value,b$value)
cohensD(c$value,d$value)
cohensD(e$value,f$value)

#Cohen'sD value
#for AC CRM= 1.63
#for VG CRM= 1.36
#For Ub CRM= 0.51




