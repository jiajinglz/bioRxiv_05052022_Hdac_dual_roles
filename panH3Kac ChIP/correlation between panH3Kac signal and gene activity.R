########################################################################
##########Differential H3Kac binding and gene expression################
########################################################################
ac_only <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/st105_panh3kac_AC_ONLY_overlapIDR_peaks.annotate", header=FALSE)
vg_only <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/st105_panh3kac_VG_ONLY_overlapIDR_peaks.annotate", header=FALSE)
ubiqutious <- read.delim("~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/st105_panh3kac_AC_VG_shared_overlapIDR_peaks.annotate", header=FALSE)

###check signal distribution###
hist(ac_only$V7)
hist(vg_only$V7)

##########Drop low signal peaks################
ac_only <- subset(ac_only, ac_only$V7> 6)
vg_only <- subset(vg_only, vg_only$V7> 6)

ac <- unique(ac_only$V14)
vg <- unique(vg_only$V14)
ubiq <- unique(ubiqutious$V7)

mydata <- list(
  AC= ac,
  VG= vg)
library("ggVennDiagram")
library(ggplot2)
p <- ggVennDiagram(mydata, label_alpha = 0, label = "count")+
  scale_fill_gradient(low = "white", high = "white")+
  scale_color_manual(values=c("tomato4", "tomato4", "tomato4", "tomato4"))
p

#ggsave("~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/st105_panh3kac_AC_VG_gene_overlap.tiff", p, dpi = 200)

#################Leave out genes in between##################

venn_info <- process_region_data(Venn(mydata))$item
ac_gene <- c(venn_info[[1]])
vg_gene <- c(venn_info[[2]])

#################Dissection data####################
tpm_Dstable <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 rna seq/TPM_table_ira_dissected_st105_tissues.txt", header= FALSE)
xtconvert <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/XENTR_10.0_GeneID_GeneName_converter.txt", header=FALSE)
geneName <- xtconvert[match(tpm_Dstable$V1, xtconvert$V2),4]
tpm_Dstable1 <- cbind(geneName, tpm_Dstable[,2:11])
colnames(tpm_Dstable1) <- c("gene_name", "ac1", "ac2", "dm1", "dm2", "lm1", "lm1", "vm1", "vm1", "vg1", "vg2")

##average Reps##
tpm_Dstable2 <- data.frame(matrix(nrow= nrow(tpm_Dstable1), ncol=4))
tpm_Dstable2[,1] <- tpm_Dstable1[,1]
tpm_Dstable2[,2] <- rowMeans(tpm_Dstable1[,2:3])
tpm_Dstable2[,3] <- rowMeans(tpm_Dstable1[,4:9])
tpm_Dstable2[,4] <- rowMeans(tpm_Dstable1[,10:11])
colnames(tpm_Dstable2) <- c("gene_name", "ac", "mz", "vg")

######only look at genes> 1 TPM########
library(rje)
tpm_Dstable3 <- subset(tpm_Dstable2, as.numeric(rowMaxs(tpm_Dstable2[, 2:4]))> 1)

########log transform###########
log_Dstable <- cbind(tpm_Dstable3$gene_name, log(tpm_Dstable3[, 2:4]+1,2))
colnames(log_Dstable) <- c("gene_name", "ac", "mz", "vg")
rownames(log_Dstable) <- log_Dstable[, 1]
log_Dstable <- log_Dstable[, -1]

########################H3Kac in AC###########################

ac_logTable <- log_Dstable[match(ac_gene, rownames(log_Dstable)),]
ac_logTable <- ac_logTable[complete.cases(ac_logTable),]
write.table(ac_logTable, file = "~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/AC_genes_germ_layer_logTPM.txt", quote = F, sep = "\t")
ac_vpTable <-stack(ac_logTable)

library(ggplot2)
p <- ggplot(ac_vpTable, aes(x= ind, y= values, fill= ind))+ 
  geom_boxplot(color="black", notch = TRUE, outlier.shape = 1)+
  stat_summary(fun=mean, geom="point", size=2, color= "red")+
  scale_fill_manual(values = c("white", "white", "white"))+
  theme_minimal()+
  theme(aspect.ratio = 1)+
  theme(axis.text = element_text(color="black"), text = element_text(size=20))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

t.test(ac_logTable$ac, ac_logTable$vg, alternative = c("two.sided"))

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/ac genes expressed in 3 germ layers.tiff', p, device = "tiff", dpi = 200)


########################H3Kac in VG###########################

vg_logTable <- log_Dstable[match(vg_gene, rownames(log_Dstable)),]
vg_logTable <- vg_logTable[complete.cases(vg_logTable),]
write.table(vg_logTable, file = "~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/VG_genes_germ_layer_logTPM.txt", quote = F, sep = "\t")
vg_vpTable <-stack(vg_logTable)

library(ggplot2)
p <- ggplot(vg_vpTable, aes(x= ind, y= values, fill= ind))+ 
  geom_boxplot(color="black", notch = TRUE, outlier.shape = 1)+
  stat_summary(fun=mean, geom="point", size=2, color= "red")+
  scale_fill_manual(values = c("white", "white", "white"))+
  theme_minimal()+
  theme(aspect.ratio = 1)+
  theme(axis.text = element_text(color="black"), text = element_text(size=20))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

t.test(vg_logTable$ac, vg_logTable$vg)

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/vg genes expressed in 3 germ layers.tiff', p, device = "tiff", dpi = 200)

################################################################################
#####################Check shared gene expression level#########################
################################################################################
ac_gene <- unique(ac_only$V14)
vg_gene <- unique(vg_only$V14)
share_gene <- unique(ubiqutious$V7)

mydata <- list(
  OP= share_gene,
  AC= ac_gene,
  VG= vg_gene)

library("ggVennDiagram")
library(ggplot2)
p <- ggVennDiagram(mydata, label_alpha = 0, label = "count")+
  scale_fill_gradient(low = "white", high = "white")

p

#ggsave("~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/st105_panh3kac_AC_VG_share_gene_overlap.tiff", p, dpi = 200)

##remove overlap genes##
venn_info <- process_region_data(Venn(mydata))$item
ac_gene <- venn_info[[2]]
vg_gene <- venn_info[[3]]
share_gene <- venn_info[[1]]

##check their regional expression##
tpm_Dstable <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 rna seq/TPM_table_ira_dissected_st105_tissues.txt", header= FALSE)
xtconvert <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/XENTR_10.0_GeneID_GeneName_converter.txt", header=FALSE)
geneName <- xtconvert[match(tpm_Dstable$V1, xtconvert$V2),4]
tpm_Dstable1 <- cbind(geneName, tpm_Dstable[,2:11])
colnames(tpm_Dstable1) <- c("gene_name", "ac1", "ac2", "dm1", "dm2", "lm1", "lm1", "vm1", "vm1", "vg1", "vg2")

##average Reps##
tpm_Dstable2 <- data.frame(matrix(nrow= nrow(tpm_Dstable1), ncol=4))
tpm_Dstable2[,1] <- tpm_Dstable1[,1]
tpm_Dstable2[,2] <- rowMeans(tpm_Dstable1[,2:3])
tpm_Dstable2[,3] <- rowMeans(tpm_Dstable1[,4:9])
tpm_Dstable2[,4] <- rowMeans(tpm_Dstable1[,10:11])
colnames(tpm_Dstable2) <- c("gene_name", "ac", "mz", "vg")


ac_gene_table <- tpm_Dstable2[match(ac_gene, tpm_Dstable2$gene_name), ]
mean(ac_gene_table$ac)
mean(ac_gene_table$mz)
mean(ac_gene_table$vg)

share_gene_table <- tpm_Dstable2[match(share_gene, tpm_Dstable2$gene_name), ]
mean(share_gene_table$ac)
mean(share_gene_table$mz)
mean(share_gene_table$vg)

vg_gene_table <- tpm_Dstable2[match(vg_gene, tpm_Dstable2$gene_name), ]
mean(vg_gene_table$ac)
mean(vg_gene_table$mz)
mean(vg_gene_table$vg)


#######Use time-course data########
gene_expression <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/time course/hdac1 regions h3k27me3 and h3k27ac/tpm_table_RdRNAseq_0_to_23hpf.txt", header=TRUE)
a <- c(1, 9)
st105_expression <- gene_expression[, a]
xtconvert <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 chip seq/time course/rna expression/XENTR_10.0_GeneID_GeneName_converter.txt", header=FALSE)
geneName <- xtconvert[match(st105_expression$gene_id, xtconvert$V2),4]
st105_expression1 <- data.frame(cbind(geneName, st105_expression[, 2]))
colnames(st105_expression1) <- c("geneName", "tpm")

ac_gene_table <- st105_expression1[match(ac_gene, st105_expression1$geneName), ]
mean(as.numeric(ac_gene_table$tpm))

share_gene_table <- st105_expression1[match(share_gene, st105_expression1$geneName), ]
mean(as.numeric(share_gene_table$tpm))

vg_gene_table <- st105_expression1[match(vg_gene, st105_expression1$geneName), ]
mean(as.numeric(vg_gene_table$tpm))

#write.table(share_gene_table$geneName, "~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/st105_panh3kac_share_genes.txt", quote = FALSE, col.names = F, row.names = F)


library(dplyr)

ac_gene_table <- ac_gene_table %>% mutate(type=c("AC_genes"))
share_gene_table <- share_gene_table %>% mutate(type=c("share_genes"))
vg_gene_table <- vg_gene_table %>% mutate(type=c("VG_genes"))

vp_table <- rbind(ac_gene_table, share_gene_table, vg_gene_table)
vp_table$tpm <- as.numeric(vp_table$tpm)

library(ggplot2)
p <- ggplot(vp_table, aes(x= type, y= log(tpm+1, 2) , fill= type))+ 
  geom_boxplot(color="black", notch = TRUE, outlier.shape = 1, width= 0.4)+
  stat_summary(fun=mean, geom="point", size=2, color= "red")+
  scale_fill_manual(values = c("white", "white", "white"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(aspect.ratio = 1)
p

ggsave('~/GDrive files/foxh1 and hdac1/panH3Kac/differential peak and gene expression/st105_panh3kac_gene_expression.tiff', p, device = "tiff", dpi = 200)

x <- subset(vp_table, vp_table$type == "AC_genes")
y <- subset(vp_table, vp_table$type == "VG_genes")
z <- subset(vp_table, vp_table$type == "share_genes")

t.test(x$tpm, z$tpm)











