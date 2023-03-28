############################################################
###############VPA AC RNAseq EdgeR DE#######################
############################################################
library(edgeR)
counts <- read.delim("~/GDrive files/foxh1 and hdac1/vpa rna/VG/st105_ut_vpa_vg_counts.txt", header=T)
rownames(counts) <- counts[,1]
counts <- counts[,-1]
colnames(counts) <- c("ctl1", "ctl2", "ctl3", "vpa1", "vpa2", "vpa3")
group <- factor(substring(colnames(counts), 1, 3))
group <- relevel (group, ref= "ctl")
rep <- factor(substring(colnames(counts), 4, 4))
y <- DGEList(counts=counts,group=group)
keep <- filterByExpr(y)
table(keep)
y<- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
design <- model.matrix(~group)
y <- estimateDisp(y, design)
plotMDS(y) ##plot PCA for all samples###

################################################################################
################################AC DE genes#####################################
################################################################################
counts <- read.delim("~/GDrive files/foxh1 and hdac1/vpa rna/VG/st105_ut_vpa_vg_counts.txt", header=T)
rownames(counts) <- counts[,1]
counts <- counts[,-1]
rep1_2 <- c(1, 2, 4, 5)
rep1_3 <- c(1, 3, 4, 6)
rep2_3 <- c(2, 3, 5, 6)
counts <- counts[, rep1_3]
#####################################################
colnames(counts) <- c("ctl1", "ctl2", "vpa1", "vpa2")
group <- factor(substring(colnames(counts), 1, 3))
group <- relevel (group, ref= "ctl")
rep <- factor(substring(colnames(counts), 4, 4))
y <- DGEList(counts=counts,group=group)
keep <- filterByExpr(y)
table(keep)
y<- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
design <- model.matrix(~group)
y <- estimateDisp(y, design)
plotMDS(y) ##plot PCA##
design <- model.matrix(~rep+rep:group)
logFC <- predFC(y, design, prior.count=1, dispersion=0.05)
cor(logFC)
design <- model.matrix(~rep+group)
rownames(design) <- colnames(y)
design
y <- estimateDisp(y, design, robust= TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmQLFit(y, design, robust = TRUE)
plotQLDisp(fit)
qlf <- glmQLFTest(fit)
topTags(qlf)
summary(decideTests(qlf))######check how many are less than 0.05 FDR####
FDR <- as.data.frame(p.adjust(qlf$table$PValue, method= "BH"))
colnames(FDR) <- c("FDR")
sum(FDR <0.05)##number should match summary above####
DElist <- as.data.frame(qlf$table)
DElist <- cbind(DElist, FDR)
DElist1 <- DElist[intersect(rownames(DElist[which(abs(DElist$logFC)> 1),]),                                               
                            rownames(DElist[which(DElist$FDR< 0.05),])),]
write.table(DElist1,file="~/GDrive files/foxh1 and hdac1/vpa rna/VG/DEgene_stat_VPA.txt", quote = FALSE, sep = "\t")
######export table######
downgenes <- DElist1[which(DElist1$logFC<0),]
write.table(downgenes,file="~/GDrive files/foxh1 and hdac1/vpa rna/VG/vg_downgenes_1_3.txt", quote = FALSE, sep = "\t")
upgenes <- DElist1[which(DElist1$logFC>0),]
write.table(upgenes,file="~/GDrive files/foxh1 and hdac1/vpa rna/VG/vg_upgenes_1_3.txt", quote = FALSE, sep = "\t")

################################################################################
#############################Compare with TSA###################################
################################################################################
ac_up_1_2 <- read.delim("~/GDrive files/foxh1 and hdac1/vpa rna/VG/vg_upgenes_1_2.txt", header=T)
ac_up_1_3 <- read.delim("~/GDrive files/foxh1 and hdac1/vpa rna/VG/vg_upgenes_1_3.txt", header=T)
ac_up_2_3 <- read.delim("~/GDrive files/foxh1 and hdac1/vpa rna/VG/vg_upgenes_2_3.txt", header=T)
ac_tsa_up <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 rna seq/endoderm/VG_upDEgene_names.txt", header=T)

ac_down_1_2 <- read.delim("~/GDrive files/foxh1 and hdac1/vpa rna/VG/vg_downgenes_1_2.txt", header=T)
ac_down_1_3 <- read.delim("~/GDrive files/foxh1 and hdac1/vpa rna/VG/vg_downgenes_1_3.txt", header=T)
ac_down_2_3 <- read.delim("~/GDrive files/foxh1 and hdac1/vpa rna/VG/vg_downgenes_2_3.txt", header=T)
ac_tsa_down <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 rna seq/endoderm/VG_downDEgene_names.txt", header=T)



summary(intersect(rownames(ac_up_1_2), rownames(ac_tsa_up)))
summary(intersect(rownames(ac_up_1_3), rownames(ac_tsa_up)))
summary(intersect(rownames(ac_up_2_3), rownames(ac_tsa_up)))

summary(intersect(rownames(ac_down_1_2), rownames(ac_tsa_down)))
summary(intersect(rownames(ac_down_1_3), rownames(ac_tsa_down)))
summary(intersect(rownames(ac_down_2_3), rownames(ac_tsa_down)))







mydata <- list(
  A= rownames(ac_down_2_3),
  B= rownames(ac_tsa_down))

library("ggVennDiagram")
library(ggplot2)
p <- ggVennDiagram(mydata, label_alpha = 0, label = "count")+
  scale_fill_gradient(low="white",high = "white")
p

vendata <- process_region_data(Venn(mydata))
share <- vendata$item[[3]]

