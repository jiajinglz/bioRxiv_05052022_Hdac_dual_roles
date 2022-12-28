################################################################################
#################Differential peak analysis plot from Homer#####################
################################################################################
####load libraries###
library(ggplot2)
library(ggrepel)
library(dplyr)
#####################
#########st9#########
st9_dfTable <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/differential peak analysis/st9_differential_peaks_result.txt", header=FALSE, comment.char="#")
st9_table <- data.frame(matrix(nrow = nrow(st9_dfTable), ncol= 3))
st9_table [, 1] <- st9_dfTable$V1
st9_table [, 2] <- log(st9_dfTable$V10, 2)
st9_table [, 3] <- st9_dfTable$V11
colnames(st9_table) <- c("st9_peaks", "Log2FC", "pValue")

st9_volTable <- st9_table %>%
  mutate(threshold = factor(case_when(Log2FC > 2 & pValue < 0.0001 ~ "up-regulated",
                                      Log2FC < -2 & pValue < 0.0001 ~ "down-regulated",
                                      TRUE ~ "n.s")))

vp <- ggplot(data=st9_volTable, aes(x=Log2FC, y=-log10(pValue))) + 
  geom_point(aes(color=threshold), alpha= 0.5, size= 2)+ 
  geom_vline(xintercept=c(-2, 2), color="purple4", lty= 2)+ 
  geom_hline(yintercept=4, color="purple4", lty= 2)+ 
  xlab("log2FC")+ 
  ylab("-log10(pValue)")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 1)+
  theme(axis.text = element_text(color="black"), text = element_text(size=9))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  scale_color_manual(name = "Genes",
                     values = c("up-regulated" = "firebrick3", 
                                "down-regulated" = "dodgerblue3", 
                                "n.s" = "grey"))
vp
#ggsave('~/GDrive files/foxh1 and hdac1/hdac1 chip seq/differential peak analysis/st9_dmso_tsa_differential_peaks.tiff', vp, device = "tiff", dpi = 200)

#####################
######st10.5#########
st105_dfTable <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/differential peak analysis/st105_differential_peaks_result.txt", header=FALSE, comment.char="#")
st105_table <- data.frame(matrix(nrow = nrow(st105_dfTable), ncol= 3))
st105_table [, 1] <- st105_dfTable$V1
summary(st105_dfTable$V10 < 0.25)
summary(st105_dfTable$V10 > 4)
st105_table [, 2] <- log(st105_dfTable$V10, 2)
st105_table [, 3] <- st105_dfTable$V11
colnames(st105_table) <- c("st105_peaks", "Log2FC", "pValue")
summary(abs(st105_table$Log2FC) > 2)
summary(abs(st105_table$Log2FC) > 2 & st105_table$pValue < 0.0001)
sig <- subset(st105_table, abs(st105_table$Log2FC) > 2 & st105_table$pValue < 0.0001)


st105_volTable <- st105_table %>%
  mutate(threshold = factor(case_when(Log2FC > 2 & pValue < 0.0001 ~ "up-regulated",
                                      Log2FC < -2 & pValue < 0.0001 ~ "down-regulated",
                                      TRUE ~ "n.s")))

vp <- ggplot(data=st105_volTable, aes(x=Log2FC, y=-log10(pValue))) + 
  geom_point(aes(color=threshold), alpha= 0.5, size= 2)+ 
  geom_vline(xintercept=c(-2, 2), color="purple4", lty= 2)+ 
  geom_hline(yintercept=4, color="purple4", lty= 2)+ 
  xlab("log2FC")+ 
  ylab("-log10(pValue)")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 1)+
  theme(axis.text = element_text(color="black"), text = element_text(size=9))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  scale_color_manual(name = "Genes",
                     values = c("up-regulated" = "firebrick3", 
                                "down-regulated" = "dodgerblue3", 
                                "n.s" = "grey"))
vp
ggsave('~/GDrive files/foxh1 and hdac1/hdac1 chip seq/differential peak analysis/st105_dmso_tsa_differential_peaks.tiff', vp, device = "tiff", dpi = 200)

################################################################################
############Differential peak analysis plot: DEseq2 UT vs amanitin##############
################################################################################
countsTable <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/differential peak analysis/st9_ut_amanitin_counts.txt")
countsTable <- as.data.frame(countsTable)
rownames(countsTable) <- make.names(countsTable[,1], unique= TRUE)
countsTable <- countsTable[,-1]
library(DESeq2)
conds <- factor(c("CTL","CTL","TRT","TRT"))
colData=data.frame(condition=conds)
dds <- DESeqDataSetFromMatrix(countsTable,colData,design=~condition)
#Perform differential expression analysis
dds2 <- DESeq(dds)
deseq_result <- results(dds2)
# identify differential expressed genes using p-value and fold change as criteria
foldchange <- rownames(deseq_result[which(abs(deseq_result$log2FoldChange)>= 1),])
adj_p <- rownames(deseq_result[which(deseq_result$padj<= 0.05),])
de_genes <- intersect(foldchange,adj_p)







