################################################################################
####################Hdac2 peak signals in Hdac1 IDR peaks#######################
################################################################################
###load libraries###
library(dplyr)
library(ggplot2)
library(RColorBrewer)

################################################################################
st9_lone <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac2/counts_st9_hdac2_lone_peaks.bed", header=FALSE)
st9_overlap <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac2/counts_st9_both_peaks.bed", header=FALSE)
st105_lone <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac2/counts_st105_hdac2_lone_peaks.bed", header=FALSE)
st105_overlap <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac2/counts_st105_both_peaks.bed", header=FALSE)

###signal= counts per kb###
st9_lone_signal <- data.frame(signal= st9_lone$V11 / (0.001*(st9_lone$V3- st9_lone$V2)))
mean(st9_lone_signal$signal)
st9_overlap_signal <- data.frame(signal= st9_overlap$V11 / (0.001*(st9_overlap$V3- st9_overlap$V2)))
mean(st9_overlap_signal$signal)
st105_lone_signal <- data.frame(signal= st105_lone$V11 / (0.001*(st105_lone$V3- st105_lone$V2)))
mean(st105_lone_signal$signal)
st105_overlap_signal <- data.frame(signal= st105_overlap$V11/ (0.001*(st105_overlap$V3- st105_overlap$V2)))
mean(st105_overlap_signal$signal)

st9_lone_signal <- st9_lone_signal %>% mutate(stage= c("st9"), type= c("lone"))
st9_overlap_signal <- st9_overlap_signal %>% mutate(stage= c("st9"), type= c("overlap"))
st105_lone_signal <- st105_lone_signal %>% mutate(stage= c("st105"), type= c("lone"))
st105_overlap_signal <- st105_overlap_signal %>% mutate(stage= c("st105"), type= c("overlap"))

plot_table <- rbind(st9_lone_signal, st9_overlap_signal, st105_lone_signal, st105_overlap_signal)

plot_table$stage <- factor(plot_table$stage, levels= c("st9", "st105"))
p <- ggplot(plot_table, aes(x= stage, y= log(signal,2)))+ 
  geom_violin(aes(fill= type), position=position_dodge(.9), trim=T, scale = "area", color= "white")+
  stat_summary(fun=mean, aes(group= type), position=position_dodge(.9), geom="point", size=2, color= "black")+
  scale_fill_manual(values=c("paleturquoise3", "royalblue"))+
  geom_boxplot(width=0.1, aes(fill= type), position=position_dodge(.9), color="grey", alpha=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 1)
p

ggsave('~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac2/hdac2_signal_on_Hdac1_peaks.tiff', p, device = "tiff", dpi = 200, width = 20, height = 20, units = c("cm"))


x <- subset(plot_table, plot_table$stage == "st105" & plot_table$type == "lone")
y <- subset(plot_table, plot_table$stage == "st105" & plot_table$type == "overlap")
t.test(log(x$signal,2), log(y$signal, 2))










