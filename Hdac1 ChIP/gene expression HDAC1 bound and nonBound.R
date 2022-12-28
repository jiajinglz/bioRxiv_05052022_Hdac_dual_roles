###########################################################################
######time-course RNA expression of Genes bound and non-bound by HDAC1#####
###########################################################################
###########################PolyA tail RNAseq###############################
hdac1_peaks <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/st9_st105_HDAC1_combined_peaks.annotation", header=FALSE)
hdac1_genes <- unique(hdac1_peaks$V7)
tpm_table <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/tpm_table_time_course_polARNAseq.txt", header=TRUE)
mean_tpm_table <- data.frame(matrix(nrow= nrow(tpm_table), ncol=9))
mean_tpm_table[,1] <- tpm_table$gene_id
mean_tpm_table[,2] <- rowMeans(tpm_table[,2:3], na.rm= TRUE)
mean_tpm_table[,3] <- rowMeans(tpm_table[,4:5], na.rm= TRUE)
mean_tpm_table[,4] <- rowMeans(tpm_table[,6:7], na.rm= TRUE)
mean_tpm_table[,5] <- rowMeans(tpm_table[,8:9], na.rm= TRUE)
mean_tpm_table[,6] <- rowMeans(tpm_table[,10:11], na.rm= TRUE)
mean_tpm_table[,7] <- rowMeans(tpm_table[,12:13], na.rm= TRUE)
mean_tpm_table[,8] <- rowMeans(tpm_table[,14:15], na.rm= TRUE)
mean_tpm_table[,9] <- rowMeans(tpm_table[,16:17], na.rm= TRUE)
xtconvert <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/XENTR_10.0_GeneID_GeneName_converter.txt", header=FALSE)
geneName <- xtconvert[match(mean_tpm_table$X1, xtconvert$V2),4]
mean_tpm_table[,1] <- geneName
colnames(mean_tpm_table) <- c("gene_name", "0hpf", "1hpf", "3hpf", "4hpf", "5hpf", "5.5hpf", "7hpf", "8hpf")
#############genes bound by HDAC1############
tpm_hdac1_gene <- subset(mean_tpm_table, mean_tpm_table$gene_name %in% hdac1_genes)
logTPM_hdac1_gene <- cbind(tpm_hdac1_gene$gene_name ,log(tpm_hdac1_gene[2:9]+1,2))
rownames(logTPM_hdac1_gene) <- logTPM_hdac1_gene[,1]
logTPM_hdac1_gene <- logTPM_hdac1_gene[,-1]
############remove any NA##############
logTPM_hdac1_gene <- logTPM_hdac1_gene[complete.cases(logTPM_hdac1_gene),]
############remove any infinity##########
logTPM_hdac1_gene <- logTPM_hdac1_gene[!is.infinite(rowSums(logTPM_hdac1_gene)),]
############remove row std=0#############
logTPM_hdac1_gene <- logTPM_hdac1_gene[!apply(logTPM_hdac1_gene , 1 , function(x) sd(x)==0 ), ]
#############genes not bound by HDAC1########
tpm_nonhdac1_gene <- subset(mean_tpm_table, !(mean_tpm_table$gene_name %in% hdac1_genes))
logTPM_nonhdac1_gene <- cbind(tpm_nonhdac1_gene$gene_name ,log(tpm_nonhdac1_gene[2:9]+1,2))
rownames(logTPM_nonhdac1_gene) <- logTPM_nonhdac1_gene[,1]
logTPM_nonhdac1_gene <- logTPM_nonhdac1_gene[,-1]
############remove any NA##############
logTPM_nonhdac1_gene <- logTPM_nonhdac1_gene[complete.cases(logTPM_nonhdac1_gene),]
############remove any infinity##########
logTPM_nonhdac1_gene <- logTPM_nonhdac1_gene[!is.infinite(rowSums(logTPM_nonhdac1_gene)),]
############remove row std=0#############
logTPM_nonhdac1_gene <- logTPM_nonhdac1_gene[!apply(logTPM_nonhdac1_gene , 1 , function(x) sd(x)==0 ), ]
###########################################################
#############heatmap############################
################################################
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
rep1 <- c(1, 3, 5)
rep2 <- c(2, 4, 6)
rep1av <- c(1, 9)
rep2av <- c(2, 10)
av <- c(1,2,9,10)
hp <- pheatmap(         logTPM_nonhdac1_gene, 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        clustering_method = "ward.D2", 
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

######################################################################
########################RiboZero RNAseq###############################
######################################################################
hdac1_peaks <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/st9_st105_HDAC1_combined_peaks.annotation", header=FALSE)
hdac1_genes <- unique(hdac1_peaks$V7)
tpm_table <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/tpm_table_RdRNAseq_zygotic_genes_0_to_8hpf.txt", header= FALSE)
xtconvert <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/XENTR_10.0_GeneID_GeneName_converter.txt", header=FALSE)
geneName <- xtconvert[match(tpm_table$V1, xtconvert$V2),4]
tpm_table1 <- cbind(geneName, tpm_table[,2:10])
colnames(tpm_table1) <- c("gene_name", "0hpf", "1hpf", "2hpf", "3hpf", "4hpf", "5hpf", "6hpf", "7hpf", "8hpf")

#################not separate hdac1 binding###############
tpm_table2 <- cbind(tpm_table1$gene_name, log(tpm_table1[2:10]+1,2))
rownames(tpm_table2) <- tpm_table2[,1]
tpm_table2 <- tpm_table2[,-1]
############remove any NA##############
tpm_table2 <- tpm_table2[complete.cases(tpm_table2),]
############remove any infinity##########
tpm_table2 <- tpm_table2[!is.infinite(rowSums(tpm_table2)),]
############transpose####################
tpm_table3 <- t(tpm_table2)

############remove row std=0#############
tpm_table2 <- tpm_table2[!apply(tpm_table2 , 1 , function(x) sd(x)==0 ), ]



#############genes bound by HDAC1############
tpm_hdac1_gene <- subset(tpm_table1, tpm_table1$gene_name %in% hdac1_genes)
logTPM_hdac1_gene <- cbind(tpm_hdac1_gene$gene_name ,log(tpm_hdac1_gene[2:10]+1,2))
rownames(logTPM_hdac1_gene) <- logTPM_hdac1_gene[,1]
logTPM_hdac1_gene <- logTPM_hdac1_gene[,-1]
############remove any NA##############
logTPM_hdac1_gene <- logTPM_hdac1_gene[complete.cases(logTPM_hdac1_gene),]
############remove any infinity##########
logTPM_hdac1_gene <- logTPM_hdac1_gene[!is.infinite(rowSums(logTPM_hdac1_gene)),]
############remove row std=0#############
logTPM_hdac1_gene <- logTPM_hdac1_gene[!apply(logTPM_hdac1_gene , 1 , function(x) sd(x)==0 ), ]
#############genes not bound by HDAC1########
tpm_nonhdac1_gene <- subset(tpm_table1, !(tpm_table1$gene_name %in% hdac1_genes))
logTPM_nonhdac1_gene <- cbind(tpm_nonhdac1_gene$gene_name ,log(tpm_nonhdac1_gene[2:10]+1,2))
rownames(logTPM_nonhdac1_gene) <- logTPM_nonhdac1_gene[,1]
logTPM_nonhdac1_gene <- logTPM_nonhdac1_gene[,-1]
############remove any NA##############
logTPM_nonhdac1_gene <- logTPM_nonhdac1_gene[complete.cases(logTPM_nonhdac1_gene),]
############remove any infinity##########
logTPM_nonhdac1_gene <- logTPM_nonhdac1_gene[!is.infinite(rowSums(logTPM_nonhdac1_gene)),]
############remove row std=0#############
logTPM_nonhdac1_gene <- logTPM_nonhdac1_gene[!apply(logTPM_nonhdac1_gene , 1 , function(x) sd(x)==0 ), ]
###########################################################
#############heatmap############################
################################################
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
rep1 <- c(1, 3, 5)
rep2 <- c(2, 4, 6)
rep1av <- c(1, 9)
rep2av <- c(2, 10)
av <- c(1,2,9,10)
hp <- pheatmap(         logTPM_hdac1_gene, 
                        scale = "column",
                        cutree_rows = 1, 
                        cluster_cols = F,
                        cluster_rows = T,
                        show_rownames = F,
                        clustering_method = "ward.D2", 
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
##################################################
#################Violin plot######################
##################################################
hdac1_peaks <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/st9_st105_HDAC1_combined_peaks.annotation", header=FALSE)
hdac1_genes <- unique(hdac1_peaks$V7)
tpm_table <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/tpm_table_RdRNAseq_zygotic_genes_0_to_8hpf.txt", header=FALSE)
xtconvert <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/XENTR_10.0_GeneID_GeneName_converter.txt", header=FALSE)
geneName <- xtconvert[match(tpm_table$V1, xtconvert$V2),4]
tpm_table1 <- cbind(geneName, tpm_table[,2:10])
colnames(tpm_table1) <- c("gene_name", "0hpf", "1hpf", "2hpf", "3hpf", "4hpf", "5hpf", "6hpf", "7hpf", "8hpf")
logTPM_table <- cbind(tpm_table1$gene_name ,log(tpm_table1[2:10]+1,2))
############remove any NA##############
logTPM_table <- logTPM_table[complete.cases(logTPM_table),]
############remove any infinity##########
logTPM_table <- logTPM_table[!is.infinite(rowSums(logTPM_table[,2:10])),]
############remove row std=0#############
logTPM_table <- logTPM_table[!apply(logTPM_table[2:10] , 1 , function(x) sd(x)==0 ), ]
colnames(logTPM_table) <- c("gene_name", "0hpf", "1hpf", "2hpf", "3hpf", "4hpf", "5hpf", "6hpf", "7hpf", "8hpf")
###############genes bound################
logTPM_hdac1 <- subset(logTPM_table, logTPM_table$gene_name %in% hdac1_genes)
logTPM_hdac1 <- logTPM_hdac1[complete.cases(logTPM_hdac1),]
###############genes not bound############
logTPM_nonhdac1 <- subset(logTPM_table, !(logTPM_table$gene_name %in% hdac1_genes))
logTPM_nonhdac1 <- logTPM_nonhdac1[complete.cases(logTPM_nonhdac1),]
##############make table for violin plot##############
library(dplyr)
logTPM_hdac1 <- logTPM_hdac1 %>% mutate(type= c("Bound") )
logTPM_nonhdac1 <- logTPM_nonhdac1 %>% mutate(type= c("Unbound") )
logTPM_both <- rbind(logTPM_hdac1, logTPM_nonhdac1)
a <- c(2,11)
logTPM_0hpf <- logTPM_both[,a] %>% mutate(hr= c("0hpf"))
colnames(logTPM_0hpf) <- c("tpm", "type", "hr")
b <- c(3,11)
logTPM_1hpf <- logTPM_both[,b] %>% mutate(hr= c("1hpf"))
colnames(logTPM_1hpf) <- c("tpm", "type", "hr")
c <- c(4,11)
logTPM_2hpf <- logTPM_both[,c] %>% mutate(hr= c("2hpf"))
colnames(logTPM_2hpf) <- c("tpm", "type", "hr")
d <- c(5,11)
logTPM_3hpf <- logTPM_both[,d] %>% mutate(hr= c("3hpf"))
colnames(logTPM_3hpf) <- c("tpm", "type", "hr")
e <- c(6,11)
logTPM_4hpf <- logTPM_both[,e] %>% mutate(hr= c("4hpf"))
colnames(logTPM_4hpf) <- c("tpm", "type", "hr")
f <- c(7,11)
logTPM_5hpf <- logTPM_both[,f] %>% mutate(hr= c("5hpf"))
colnames(logTPM_5hpf) <- c("tpm", "type", "hr")
g <- c(8,11)
logTPM_6hpf <- logTPM_both[,g] %>% mutate(hr= c("6hpf"))
colnames(logTPM_6hpf) <- c("tpm", "type", "hr")
h <- c(9,11)
logTPM_7hpf <- logTPM_both[,h] %>% mutate(hr= c("7hpf"))
colnames(logTPM_7hpf) <- c("tpm", "type", "hr")
i <- c(10,11)
logTPM_8hpf <- logTPM_both[,i] %>% mutate(hr= c("8hpf"))
colnames(logTPM_8hpf) <- c("tpm", "type", "hr")

logTPM_vpTable <- rbind(logTPM_0hpf, logTPM_1hpf, logTPM_2hpf, 
                        logTPM_3hpf, logTPM_4hpf, logTPM_5hpf,
                        logTPM_6hpf, logTPM_7hpf, logTPM_8hpf)

logTPM_vpTable <- rbind(logTPM_0hpf, logTPM_4hpf, logTPM_5hpf,
                        logTPM_6hpf, logTPM_7hpf, logTPM_8hpf)

library(ggplot2)
p <- ggplot(logTPM_vpTable, aes(x= hr, y= tpm, fill= type))+ 
  geom_boxplot(outlier.shape = NA)+
  ylim(-0.5, 1.5)+
  scale_fill_manual(values=c("lightskyblue", "palegoldenrod"))+
  geom_hline(yintercept = 0, linetype='dotted', col = 'red')+
  theme_minimal()+
  theme(aspect.ratio = 3/5)+  
  theme(axis.text = element_text(color="black"), text = element_text(size=9))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")
p
ggsave('G:/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/time-course RdRNA of HDAC1 binding--LEGEND.tiff', p, device = "tiff", dpi = 500)















