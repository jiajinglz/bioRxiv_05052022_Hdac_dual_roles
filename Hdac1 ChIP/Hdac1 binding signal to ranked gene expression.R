nCounts <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/normalized_counts_timeCourse.txt", header=FALSE)
avgCounts <- data.frame(matrix(nrow = nrow(nCounts), ncol = 5))
avgCounts[,1] <- nCounts[,1]
avgCounts[,2] <- round(rowMeans(nCounts[,2:3]))
avgCounts[,3] <- round(rowMeans(nCounts[,4:5]))
avgCounts[,4] <- round(rowMeans(nCounts[,6:7]))
avgCounts[,5] <- round(rowMeans(nCounts[,8:9]))
write.table(avgCounts, "C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/average_normalized_counts_timeCourse.bed",
            row.names = F, col.names = F, quote = F, sep = "\t")
#################st9################
st9_signal <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/st9_XENTR_10.0_GCF_ProteinCoding_promoter+-2kb_counts.txt", header=FALSE)
st9_geneOrder <- avgCounts[order(avgCounts$X4),]
st9_signal_geneCounts <- na.omit(st9_signal[match(st9_geneOrder$X1, st9_signal$V4),])
st9_plot_table <- st9_signal_geneCounts[,6:7]
st9_plot_table$V6 <- c(1:21819)
rownames(st9_plot_table) <- st9_signal_geneCounts$V4
#######replace 0 with 1 for log transformation###########
st9_plot_table[st9_plot_table == 0] <- 1
summary(st9_plot_table)
library(ggplot2)
# Basic dot plot
p9 <- ggplot(st9_plot_table, aes(x=V6, y=V7)) + 
      geom_point(size= 0.5, color="steelblue4", alpha =1/10)+
      scale_x_log10()+
      scale_y_log10()+
      annotation_logticks()+
      geom_smooth(method="loess", color= "tomato", se= F)+
      theme_minimal()+
      theme(aspect.ratio = 1)
p9
ggsave('C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/st9_hdac1_signal_to_ranked_geneExpression_5.5hpf.tiff', p9, device = "tiff", dpi = 500)
##################st10.5########################
st105_signal <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/st105_XENTR_10.0_GCF_ProteinCoding_promoter+-2kb_counts.txt", header=FALSE)
st105_geneOrder <- avgCounts[order(avgCounts$X5),]
st105_signal_geneCounts <- na.omit(st105_signal[match(st105_geneOrder$X1, st105_signal$V4),])
st105_plot_table <- st105_signal_geneCounts[,6:7]
st105_plot_table$V6 <- c(1:21819)
rownames(st105_plot_table) <- st105_signal_geneCounts$V4
#######replace 0 with 1 for log transformation###########
st105_plot_table[st105_plot_table == 0] <- 1
summary(st105_plot_table)
library(ggplot2)
# Basic dot plot
p105 <- ggplot(st105_plot_table, aes(x=V6, y=V7)) + 
  geom_point(size= 0.5, color="steelblue4", alpha= 1/10)+
  scale_x_log10()+
  scale_y_log10()+
  annotation_logticks()+
  geom_smooth(method="loess", color= "tomato", se= F)+
  theme_minimal()+
  theme(aspect.ratio = 1)
p105
ggsave('C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/st105_hdac1_signal_to_ranked_geneExpression_5.5hpf.tiff', p105, device = "tiff", dpi = 500)

