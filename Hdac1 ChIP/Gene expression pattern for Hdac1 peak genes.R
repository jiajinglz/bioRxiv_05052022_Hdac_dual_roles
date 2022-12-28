################################################################################
###############Hdac1 peak associated genes expression pattern###################
################################################################################

st8_unique<- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/st8_unique_IDR0.05_peaks.annotate", header=FALSE)
st9_unique<- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/st9_unique_IDR0.05_peaks.annotate", header=FALSE)
st105_unique<- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/st105_unique_IDR0.05_peaks.annotate", header=FALSE)
non_unique<- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/st8_9_105_anyTwo_IDR0.05_peaks.annotate", header=FALSE)

st8_gene <- unique(st8_unique$V14)
st9_gene <- unique(st9_unique$V14)
st105_gene <- unique(st105_unique$V14)
share_gene <- unique(non_unique$V7)

mydata <- list(
  A= st8_gene, 
  B= st9_gene, 
  C= share_gene, 
  D= st105_gene)
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library(ggplot2)
p <- ggVennDiagram(mydata, label_alpha = 0, label = "count")+
  scale_fill_gradient(low="white",high = "white")
p

#ggsave('~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/hdac1_timecourse_cluster_abc_genes_Venn.tiff', p, device = "tiff", dpi = 200)

vendata <- process_region_data(Venn(mydata))
st8_gene <- vendata$item[[1]]
st9_gene <- vendata$item[[2]]
share_gene <- vendata$item[[3]]

tpm_zygGene<- read.delim("~/GDrive files/foxh1 and hdac1/time course rna/tpm_table_RdRNAseq_0_to_23hpf_FUNCTIONAL_zygotic_genes.txt", header=TRUE)
tpm_zygGene <- tpm_zygGene[, 1:12]

#tpm_st8_gene <- tpm_allGene[match(st8_gene, tpm_allGene$gene_id),]
tpm_st8_gene <- tpm_zygGene[match(st8_gene, tpm_zygGene$gene_name),]
tpm_st8_gene <- tpm_st8_gene[complete.cases(tpm_st8_gene),]
colnames(tpm_st8_gene) <- c("gene_id", "0hr", "1hr", "2hr", "3hr", "4hr", "5hr",
                            "6hr", "7hr", "8hr", "9hr", "10hr")
#tpm_st9_gene <- tpm_allGene[match(st9_gene, tpm_allGene$gene_id),]
tpm_st9_gene <- tpm_zygGene[match(st9_gene, tpm_zygGene$gene_name),]
tpm_st9_gene <- tpm_st9_gene[complete.cases(tpm_st9_gene),]
colnames(tpm_st9_gene) <- c("gene_id", "0hr", "1hr", "2hr", "3hr", "4hr", "5hr",
                            "6hr", "7hr", "8hr", "9hr", "10hr")
#tpm_share_gene <- tpm_allGene[match(share_gene, tpm_allGene$gene_id),]
tpm_share_gene <- tpm_zygGene[match(share_gene, tpm_zygGene$gene_name),]
tpm_share_gene <- tpm_share_gene[complete.cases(tpm_share_gene),]
colnames(tpm_share_gene) <- c("gene_id", "0hr", "1hr", "2hr", "3hr", "4hr", "5hr",
                            "6hr", "7hr", "8hr", "9hr", "10hr")
#tpm_st105_gene <- tpm_allGene[match(st105_gene, tpm_allGene$gene_id),]
tpm_st105_gene <- tpm_zygGene[match(st105_gene, tpm_zygGene$gene_name),]
tpm_st105_gene <- tpm_st105_gene[complete.cases(tpm_st105_gene),]
colnames(tpm_st105_gene) <- c("gene_id", "0hr", "1hr", "2hr", "3hr", "4hr", "5hr",
                              "6hr", "7hr", "8hr", "9hr", "10hr")



tpm_st8_gene <- stack(tpm_st8_gene[, 2:12])
tpm_st9_gene <- stack(tpm_st9_gene[, 2:12])
tpm_share_gene <- stack(tpm_share_gene[, 2:12])
tpm_st105_gene <- stack(tpm_st105_gene[, 2:12])



library(dplyr)
tpm_st8_gene <- tpm_st8_gene %>% mutate(type=c("cluster a"))
tpm_st9_gene <- tpm_st9_gene %>% mutate(type=c("cluster b"))
tpm_share_gene <- tpm_share_gene %>% mutate(type=c("cluster c"))
tpm_st105_gene <- tpm_st105_gene %>% mutate(type=c("cluster d"))
tpm_table <- rbind(tpm_st8_gene, tpm_st9_gene, tpm_share_gene, tpm_st105_gene)
colnames(tpm_table) <- c("TPM", "time", "type")

library(ggplot2)
p <- ggplot(tpm_table, aes(x= time, y= log(TPM+1, 2), fill= type))+ 
  geom_boxplot(outlier.color = "grey", outlier.shape = NA, coef= 0.5)+
  scale_fill_manual(values=alpha(c("cyan", "skyblue1", "skyblue3", "steelblue"), 1))+
  geom_hline(yintercept=0,linetype=2, color= c("red"))+
  scale_x_discrete(labels= c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))+
  ylim(0, 5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 3/5)

p

ggsave('~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/hdac1_timecourse_clusters_zygotic_genes.tiff', p, device = "tiff", width= 20, height= 20, units= c("cm"), dpi = 200)

##notice the transient cluster of st8 expression##

x <- subset(tpm_table, tpm_table$time == "3hr"& tpm_table$type == "cluster c")
y <- subset(tpm_table, tpm_table$time == "4hr"& tpm_table$type == "cluster c")
z <- subset(tpm_table, tpm_table$time == "5hr"& tpm_table$type == "cluster c")

t.test(log(x$TPM+1,2), log(y$TPM+1,2))
t.test(log(y$TPM+1,2), log(z$TPM+1,2))

#####################################################
####################Heatmap##########################
#####################################################
hp_st8_table <- tpm_st8_gene[apply(tpm_st8_gene[, 5:10], 1, sd, na.rm=TRUE) != 0,]
library(pheatmap)
library(ggplot2)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
a_hp <- pheatmap(       log(hp_st8_table[, 5:10]+1, 2), 
                             scale = "row",
                             cutree_rows = 1, 
                             cluster_cols = F, 
                             show_rownames = F,
                             clustering_method = "ward.D2", 
                             border_color = NA,
                             cellwidth = 30,
                             cellheight = 0.2,
                             treeheight_row = 0,
                             col=colorRampPalette(brewer.pal(4,"Blues"))(255))

ggsave('~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/st8_gene_hp.tiff', a_hp, device = "tiff", dpi = 300)

hp_st9_table <- tpm_st9_gene[apply(tpm_st9_gene[, 5:10], 1, sd, na.rm=TRUE) != 0,]
library(pheatmap)
library(ggplot2)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
b_hp <- pheatmap(       log(hp_st9_table[, 5:10]+1, 2), 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        clustering_method = "ward.D2", 
                        border_color = NA,
                        cellwidth = 30,
                        cellheight = 0.2,
                        treeheight_row = 0,
                        col=colorRampPalette(brewer.pal(4,"Blues"))(255))

ggsave('~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/st9_gene_hp.tiff', b_hp, device = "tiff", dpi = 300)


hp_share_table <- tpm_share_gene[apply(tpm_share_gene[, 5:10], 1, sd, na.rm=TRUE) != 0,]
library(pheatmap)
library(ggplot2)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
c_hp <- pheatmap(       log(hp_share_table[, 5:10]+1, 2), 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        clustering_method = "ward.D2", 
                        border_color = NA,
                        cellwidth = 30,
                        cellheight = 0.2,
                        treeheight_row = 0,
                        col=colorRampPalette(brewer.pal(4,"Blues"))(255))

ggsave('~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/share_gene_hp.tiff', c_hp, device = "tiff", dpi = 300)


hp_st105_table <- tpm_st105_gene[apply(tpm_st105_gene[, 5:10], 1, sd, na.rm=TRUE) != 0,]
library(pheatmap)
library(ggplot2)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
d_hp <- pheatmap(       log(hp_st105_table[, 5:10]+1, 2), 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        clustering_method = "ward.D2", 
                        border_color = NA,
                        cellwidth = 30,
                        cellheight = 0.2,
                        treeheight_row = 0,
                        col=colorRampPalette(brewer.pal(4,"Blues"))(255))

ggsave('~/GDrive files/foxh1 and hdac1/hdac1 chip seq/hdac1 peaks associated genes/st105_gene_hp.tiff', d_hp, device = "tiff", dpi = 300)


