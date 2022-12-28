###################################################################
################gene annotation in each cluster####################
###################################################################
cluster1 <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/st9_st105_hdac1_h3k27ac_h3k27me3_overlap_peaks.annotate", header=FALSE)
cluster2 <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/st9_st105_hdac1_only_h3k27me3_peaks.annotate", header=FALSE)
cluster3 <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/st9_st105_hdac1_only_h3k27ac_peaks.annotate", header=FALSE)
cluster4 <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/st9_st105_hdac1_neither_h3k27ac_h3k27me3_peaks.annotate", header=FALSE)

############assignment is valid only within 10kb############
cluster1 <- subset(cluster1, abs(cluster1$V10)< 10000)
cluster2 <- subset(cluster2, abs(cluster2$V10)< 10000)
cluster3 <- subset(cluster3, abs(cluster3$V10)< 10000)
cluster4 <- subset(cluster4, abs(cluster4$V10)< 10000)


cluster1_gene <- unique(cluster1$V7)
cluster2_gene <- unique(cluster2$V7)
cluster3_gene <- unique(cluster3$V7)
cluster4_gene <- unique(cluster4$V7)


mydata <- list(
  I= cluster1_gene,
  II= cluster2_gene,
  III= cluster3_gene)
##if (!require(devtools)) install.packages("devtools")
##devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library(ggplot2)
p <- ggVennDiagram(mydata, label_alpha = 0, label = "count")+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values=c("tomato4", "tomato4", "tomato4", "tomato4"))
p

#ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/Venn_of_genes_in_3_clusters.tiff', p, device = "tiff", dpi = 200)

venn_info <- process_region_data(Venn(mydata))$item
library(RJSONIO)
exportJSON <- toJSON(venn_info)
#write(exportJSON,"/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/venn_gene_classes_info_sum.txt")

class1 <- c(venn_info[[1]], venn_info[[4]], venn_info[[5]], venn_info[[6]], venn_info[[7]])
#write.table(class1, file = "/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/class1_gene_both_h3k27.txt", quote = F, row.names = F, col.names = F)
class2 <- c(venn_info[[2]])
#write.table(class2, file = "/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/class2_gene_h3k27me3_only.txt", quote = F, row.names = F, col.names = F)
class3 <- c(venn_info[[3]])
#write.table(class3, file = "/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/class3_gene_h3k27ac_only.txt", quote = F, row.names = F, col.names = F)


###################################################################
################gene expression in each cluster####################
###################################################################
tpm_tpm_zygoGene_table <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/time course rna/tpm_table_RdRNAseq_0_to_23hpf_FUNCTIONAL_zygotic_genes.txt", header= T)
tpm_Rdtable1 <- tpm_zygoGene_table
colnames(tpm_Rdtable1) <- c("gene_name", "0hpf", "1hpf", "2hpf", "3hpf", "4hpf", "5hpf", "6hpf", 
                            "7hpf", "8hpf", "9hpf", "10hpf", "11hpf", "12hpf", "13hpf", "14hpf", "15hpf",
                            "16hpf", "17hpf", "18hpf", "19hpf", "20hpf", "21hpf", "22hpf", "23hpf")
logTPM_Rdtable <- cbind(tpm_Rdtable1$gene_name ,log(tpm_Rdtable1[2:25]+1,2))
colnames(logTPM_Rdtable) <- c("gene_name", "0hpf", "1hpf", "2hpf", "3hpf", "4hpf", "5hpf", "6hpf", 
                            "7hpf", "8hpf", "9hpf", "10hpf", "11hpf", "12hpf", "13hpf", "14hpf", "15hpf",
                            "16hpf", "17hpf", "18hpf", "19hpf", "20hpf", "21hpf", "22hpf", "23hpf")

logTPM_class1 <- logTPM_Rdtable[match(class1, logTPM_Rdtable$gene_name),]
logTPM_class1<- logTPM_class1[complete.cases(logTPM_class1 ),]
logTPM_class2 <- logTPM_Rdtable[match(class2, logTPM_Rdtable$gene_name),]
logTPM_class2<- logTPM_class2[complete.cases(logTPM_class2 ),]
logTPM_class3 <- logTPM_Rdtable[match(class3, logTPM_Rdtable$gene_name),]
logTPM_class3<- logTPM_class3[complete.cases(logTPM_class3 ),]


################################################################################
################three class gene expression at gastrulation#####################
#########################Functional zygotic genes only##########################
library(dplyr)
logTPM_class1 <- logTPM_class1 %>% mutate(type= c("ClassI") )
logTPM_class2 <- logTPM_class2 %>% mutate(type= c("ClassII") )
logTPM_class3 <- logTPM_class3 %>% mutate(type= c("ClassIII") )
logTPM_all <- rbind(logTPM_class1, logTPM_class2, logTPM_class3)

a <- c(2,26)
logTPM_0hpf <- logTPM_all[,a] %>% mutate(hr= c("0hpf"))
colnames(logTPM_0hpf) <- c("tpm", "type", "hr")
b <- c(5,26)
logTPM_3hpf <- logTPM_all[,b] %>% mutate(hr= c("3hpf"))
colnames(logTPM_3hpf) <- c("tpm", "type", "hr")
c <- c(6,26)
logTPM_4hpf <- logTPM_all[,c] %>% mutate(hr= c("4hpf"))
colnames(logTPM_4hpf) <- c("tpm", "type", "hr")
d <- c(7,26)
logTPM_5hpf <- logTPM_all[,d] %>% mutate(hr= c("5hpf"))
colnames(logTPM_5hpf) <- c("tpm", "type", "hr")
e <- c(8,26)
logTPM_6hpf <- logTPM_all[,e] %>% mutate(hr= c("6hpf"))
colnames(logTPM_6hpf) <- c("tpm", "type", "hr")
f <- c(9,26)
logTPM_7hpf <- logTPM_all[,f] %>% mutate(hr= c("7hpf"))
colnames(logTPM_7hpf) <- c("tpm", "type", "hr")
g <- c(10,26)
logTPM_8hpf <- logTPM_all[,g] %>% mutate(hr= c("8hpf"))
colnames(logTPM_8hpf) <- c("tpm", "type", "hr")
h <- c(11,26)
logTPM_9hpf <- logTPM_all[,h] %>% mutate(hr= c("9hpf"))
colnames(logTPM_9hpf) <- c("tpm", "type", "hr")
i <- c(12,26)
logTPM_10hpf <- logTPM_all[,i] %>% mutate(hr= c("10hpf"))
colnames(logTPM_10hpf) <- c("tpm", "type", "hr")


logTPM_vpTable <- rbind(logTPM_0hpf, logTPM_3hpf, logTPM_4hpf, 
                        logTPM_5hpf, logTPM_6hpf, logTPM_7hpf,
                        logTPM_8hpf, logTPM_9hpf, logTPM_10hpf)

x_order <- c("0hpf", "3hpf", "4hpf", "5hpf", "6hpf", "7hpf", "8hpf", "9hpf", "10hpf")

library(ggplot2)
p <- ggplot(logTPM_vpTable, aes(x= factor(hr, level=x_order), y= tpm, fill= type))+ 
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("palegoldenrod", "tomato", "darkseagreen2"))+
  geom_hline(yintercept = 0, linetype='dotted', col = 'red')+
  scale_x_discrete(labels= c("0", "3", "4", "5", "6", "7", "8", "9", "10"))+
  ylim(0, 6)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 3/5)

p

ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/upto gastrula zygotic genes on different gene classes.tiff', p, device = "tiff", dpi = 200)

#######################################################################
############boxplot trendline over 23 hpf for each cluster#############
########################Functional Zygotic genes#######################
tpm_zygoGene_table <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/time course rna/tpm_table_RdRNAseq_0_to_23hpf_FUNCTIONAL_zygotic_genes.txt", header= T)
tpm_Rdtable1 <- tpm_zygoGene_table
colnames(tpm_Rdtable1) <- c("gene_name", "0hpf", "1hpf", "2hpf", "3hpf", "4hpf", "5hpf", "6hpf", 
                            "7hpf", "8hpf", "9hpf", "10hpf", "11hpf", "12hpf", "13hpf", "14hpf", "15hpf",
                            "16hpf", "17hpf", "18hpf", "19hpf", "20hpf", "21hpf", "22hpf", "23hpf")
logTPM_Rdtable <- cbind(tpm_Rdtable1$gene_name ,log(tpm_Rdtable1[2:25]+1,2))
colnames(logTPM_Rdtable) <- c("gene_name", "0hpf", "1hpf", "2hpf", "3hpf", "4hpf", "5hpf", "6hpf", 
                              "7hpf", "8hpf", "9hpf", "10hpf", "11hpf", "12hpf", "13hpf", "14hpf", "15hpf",
                              "16hpf", "17hpf", "18hpf", "19hpf", "20hpf", "21hpf", "22hpf", "23hpf")

logTPM_class1 <- logTPM_Rdtable[match(class1, logTPM_Rdtable$gene_name),]
logTPM_class1<- logTPM_class1[complete.cases(logTPM_class1 ),]
logTPM_class2 <- logTPM_Rdtable[match(class2, logTPM_Rdtable$gene_name),]
logTPM_class2<- logTPM_class2[complete.cases(logTPM_class2 ),]
logTPM_class3 <- logTPM_Rdtable[match(class3, logTPM_Rdtable$gene_name),]
logTPM_class3<- logTPM_class3[complete.cases(logTPM_class3 ),]

library(dplyr)
sc <- c(2, 5, 6, 7, 8, 9, 10, 12, 16, 19, 22, 25)
logTPM_class1 <- stack(logTPM_class1[, sc])
logTPM_class1 <- logTPM_class1 %>% mutate(type= c("ClassI") )
logTPM_class2 <- stack(logTPM_class2[, sc])
logTPM_class2 <- logTPM_class2 %>% mutate(type= c("ClassII") )
logTPM_class3 <- stack(logTPM_class3[, sc])
logTPM_class3 <- logTPM_class3 %>% mutate(type= c("ClassIII") )
logTPM_all_class <- rbind(logTPM_class1, logTPM_class2, logTPM_class3)

library(ggplot2)
p <- ggplot(logTPM_all_class, aes(x= ind, y= values, fill= type))+ 
  geom_boxplot(outlier.shape = NA, coef= 0)+
  scale_fill_manual(values=alpha(c("palegoldenrod", "tomato", "darkseagreen2"), 0.6))+
  geom_hline(yintercept = 0, linetype='dotted', col = 'red')+
  geom_vline(xintercept = 4, linetype='dotted', col = 'black')+
  geom_vline(xintercept = 6, linetype='dotted', col = 'black')+
  ylim(0, 6)+
  scale_x_discrete(labels= c("0", "3", "4", "5", "6", "7", "8", "10", "14", "17", "20", "23"))+
  geom_smooth(data= subset(logTPM_all_class, type == "ClassI"), se=FALSE, aes(group= 1), color= "palegoldenrod")+
  geom_smooth(data= subset(logTPM_all_class, type == "ClassII"), se=FALSE, aes(group= 1), color= "tomato")+
  geom_smooth(data= subset(logTPM_all_class, type == "ClassIII"), se=FALSE, aes(group= 1), color= "darkseagreen2")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 3/5)

p
ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/time-course_Rd_tpm_TrendLine_zygotic_genes.tiff', p, device = "tiff", dpi = 200)


#######################################################################
############boxplot trendline over 23 hpf for each cluster#############
###############################ALL genes###############################
tpm_allGene_table <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/time course rna/tpm_table_RdRNAseq_0_to_23hpf.txt", header= T)
library(stringr)
gene_name <- tpm_allGene_table$gene_id %>% str_replace("gene-", "")
tpm_allGene_table[ , 1] <- gene_name
tpm_Rdtable1 <- tpm_allGene_table
colnames(tpm_Rdtable1) <- c("gene_name", "0hpf", "1hpf", "2hpf", "3hpf", "4hpf", "5hpf", "6hpf", 
                            "7hpf", "8hpf", "9hpf", "10hpf", "11hpf", "12hpf", "13hpf", "14hpf", "15hpf",
                            "16hpf", "17hpf", "18hpf", "19hpf", "20hpf", "21hpf", "22hpf", "23hpf")
logTPM_Rdtable <- cbind(tpm_Rdtable1$gene_name ,log(tpm_Rdtable1[2:25]+1,2))
colnames(logTPM_Rdtable) <- c("gene_name", "0hpf", "1hpf", "2hpf", "3hpf", "4hpf", "5hpf", "6hpf", 
                              "7hpf", "8hpf", "9hpf", "10hpf", "11hpf", "12hpf", "13hpf", "14hpf", "15hpf",
                              "16hpf", "17hpf", "18hpf", "19hpf", "20hpf", "21hpf", "22hpf", "23hpf")

logTPM_class1 <- logTPM_Rdtable[match(class1, logTPM_Rdtable$gene_name),]
logTPM_class1<- logTPM_class1[complete.cases(logTPM_class1 ),]
logTPM_class2 <- logTPM_Rdtable[match(class2, logTPM_Rdtable$gene_name),]
logTPM_class2<- logTPM_class2[complete.cases(logTPM_class2 ),]
logTPM_class3 <- logTPM_Rdtable[match(class3, logTPM_Rdtable$gene_name),]
logTPM_class3<- logTPM_class3[complete.cases(logTPM_class3 ),]

library(dplyr)
sc <- c(2, 5, 6, 7, 8, 9, 10, 12, 16, 19, 22, 25)
logTPM_class1 <- stack(logTPM_class1[, sc])
logTPM_class1 <- logTPM_class1 %>% mutate(type= c("ClassI") )
logTPM_class2 <- stack(logTPM_class2[, sc])
logTPM_class2 <- logTPM_class2 %>% mutate(type= c("ClassII") )
logTPM_class3 <- stack(logTPM_class3[, sc])
logTPM_class3 <- logTPM_class3 %>% mutate(type= c("ClassIII") )
logTPM_all_class <- rbind(logTPM_class1, logTPM_class2, logTPM_class3)

library(ggplot2)
p <- ggplot(logTPM_all_class, aes(x= ind, y= values, fill= type))+ 
  geom_boxplot(outlier.shape = NA, coef= 0)+
  scale_fill_manual(values=alpha(c("palegoldenrod", "tomato", "darkseagreen2"), 0.6))+
  geom_hline(yintercept = 0, linetype='dotted', col = 'red')+
  geom_vline(xintercept = 4, linetype='dotted', col = 'black')+
  geom_vline(xintercept = 6, linetype='dotted', col = 'black')+
  scale_x_discrete(labels= c("0", "3", "4", "5", "6", "7", "8", "10", "14", "17", "20", "23"))+
  ylim(0, 8)+
  geom_smooth(data= subset(logTPM_all_class, type == "ClassI"), se=FALSE, aes(group= 1), color= "palegoldenrod")+
  geom_smooth(data= subset(logTPM_all_class, type == "ClassII"), se=FALSE, aes(group= 1), color= "tomato")+
  geom_smooth(data= subset(logTPM_all_class, type == "ClassIII"), se=FALSE, aes(group= 1), color= "darkseagreen2")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 3/5)

p

ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/time-course_Rd_tpm_TrendLine_all_genes.tiff', p, device = "tiff", dpi = 200)


###################################################################
##############dissected st10.5 genes in each cluster###############
###################################################################
tpm_Dstable <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/TPM_table_ira_dissected_st105_tissues.txt", header= FALSE)
library(stringr)
gene_name <- tpm_Dstable$V1 %>% str_replace("gene-", "")
tpm_Dstable[ , 1] <- gene_name
tpm_Dstable1 <- tpm_Dstable
colnames(tpm_Dstable1) <- c("gene_name", "ac1", "ac2", "dm1", "dm2", "lm1", "lm1", "vm1", "vm1", "vg1", "vg2")
##average Reps##
tpm_Dstable2 <- data.frame(matrix(nrow= nrow(tpm_Dstable1), ncol=6))
tpm_Dstable2[,1] <- tpm_Dstable1[,1]
tpm_Dstable2[,2] <- rowMeans(tpm_Dstable1[,2:3])
tpm_Dstable2[,3] <- rowMeans(tpm_Dstable1[,4:5])
tpm_Dstable2[,4] <- rowMeans(tpm_Dstable1[,6:7])
tpm_Dstable2[,5] <- rowMeans(tpm_Dstable1[,8:9])
tpm_Dstable2[,6] <- rowMeans(tpm_Dstable1[,10:11])
colnames(tpm_Dstable2) <- c("gene_name", "ac", "dm", "lm", "vm", "vg")
tpmDs_class1 <- tpm_Dstable2[match(class1, tpm_Dstable2$gene_name),]
tpmDs_class1 <- tpmDs_class1[complete.cases(tpmDs_class1),]
tpmDs_class1 <- cbind(tpmDs_class1$gene_name,log(tpmDs_class1[, 2:6]+1,2))
colnames(tpmDs_class1) <- c("gene_name", "ac", "dm", "lm", "vm", "vg")
tpmDs_class2 <- tpm_Dstable2[match(class2, tpm_Dstable2$gene_name),]
tpmDs_class2 <- tpmDs_class2[complete.cases(tpmDs_class2),]
tpmDs_class2 <- cbind(tpmDs_class2$gene_name,log(tpmDs_class2[, 2:6]+1,2))
colnames(tpmDs_class2) <- c("gene_name", "ac", "dm", "lm", "vm", "vg")
tpmDs_class3 <- tpm_Dstable2[match(class3, tpm_Dstable2$gene_name),]
tpmDs_class3 <- tpmDs_class3[complete.cases(tpmDs_class3),]
tpmDs_class3 <- cbind(tpmDs_class3$gene_name,log(tpmDs_class3[, 2:6]+1,2))
colnames(tpmDs_class3) <- c("gene_name", "ac", "dm", "lm", "vm", "vg")

############################################
##############find variance#################
############################################
############################################
sd_class1 <- data.frame(apply(tpmDs_class1[, 2:6], 1, function(x) sd(x)))
sd_class1_table <- cbind(tpmDs_class1$gene_name, sd_class1)
colnames(sd_class1_table) <- c("gene_name", "std")
sd_class2 <- data.frame(apply(tpmDs_class2[, 2:6], 1, function(x) sd(x)))
sd_class2_table <- cbind(tpmDs_class2$gene_name, sd_class2)
colnames(sd_class2_table) <- c("gene_name", "std")
sd_class3 <- data.frame(apply(tpmDs_class3[, 2:6], 1, function(x) sd(x)))
sd_class3_table <- cbind(tpmDs_class3$gene_name, sd_class3)
colnames(sd_class3_table) <- c("gene_name", "std")

library(dplyr)
sd_class1_table <- sd_class1_table %>% mutate(type= c("ClassI") )
sd_class2_table <- sd_class2_table %>% mutate(type= c("ClassII") )
sd_class3_table <- sd_class3_table %>% mutate(type= c("ClassIII") )
sd_allClass_table <- rbind(sd_class1_table, sd_class2_table, sd_class3_table)
sd_allClass_table <- sd_allClass_table[complete.cases(sd_allClass_table),]
#####################################
############Violin Plot##############
#####################################
library(ggplot2)
p <- ggplot(sd_allClass_table[, 2:3], aes(x= type, y= std, fill= type))+ 
  geom_violin(trim=T, scale = "area", color= "white")+
  stat_summary(fun=mean, geom="point", size=1, color= "red")+
  scale_fill_manual(values=alpha(c("palegoldenrod", "tomato", "darkseagreen2"), 0.7))+
  geom_boxplot(width=0.1, color="grey", alpha=0.5, outlier.shape = 1, outlier.size = 0.5)+
  geom_hline(yintercept = 0, linetype='dotted', col = 'red')+
  ylim(0, 2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 1)

p

ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/std of dissection all-gene TPM on different clusters.tiff', p, device = "tiff", dpi = 200)

t.test(sd_class1_table$std, sd_class2_table$std, alternative = "greater")
t.test(sd_class1_table$std, sd_class3_table$std, alternative = "greater")
t.test(sd_class2_table$std, sd_class3_table$std)



