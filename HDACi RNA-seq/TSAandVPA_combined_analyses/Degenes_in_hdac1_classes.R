################################################################################
##############Spatial distribution of DEgenes by TSA and VPA####################
################################################################################
#LOAD PACKAGES#
library(stringr)
library(ggplot2)
library(pheatmap)
library("RColorBrewer")
############################Import DEgenes######################################
tsa_ac_upgene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/TSA_AC_upDEgene_names.txt", header=F)
tsa_ac_downgene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/TSA_AC_downDEgene_names.txt", header=F)
tsa_vg_upgene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/TSA_VG_upDEgene_names.txt", header=F)
tsa_vg_downgene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/TSA_VG_downDEgene_names.txt", header=F)

vpa_ac_upgene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/VPA_AC_upDEgenes.txt", header=T)
vpa_ac_upgene <- rownames(vpa_ac_upgene)
vpa_ac_upgene <- str_remove(vpa_ac_upgene, "gene-")
vpa_ac_upgene <- data.frame(V1= vpa_ac_upgene)

vpa_ac_downgene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/VPA_AC_downDEgenes.txt", header=T)
vpa_ac_downgene <- rownames(vpa_ac_downgene)
vpa_ac_downgene <- str_remove(vpa_ac_downgene, "gene-")
vpa_ac_downgene <- data.frame(V1= vpa_ac_downgene)

vpa_vg_upgene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/VPA_VG_upDEgenes.txt", header=T)
vpa_vg_upgene <- rownames(vpa_vg_upgene)
vpa_vg_upgene <- str_remove(vpa_vg_upgene, "gene-")
vpa_vg_upgene <- data.frame(V1= vpa_vg_upgene)

vpa_vg_downgene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/VPA_VG_downDEgenes.txt", header=T)
vpa_vg_downgene <- rownames(vpa_vg_downgene)
vpa_vg_downgene <- str_remove(vpa_vg_downgene, "gene-")
vpa_vg_downgene <- data.frame(V1= vpa_vg_downgene)

####find co-regulated genes####
co_ac_upgene <- intersect(tsa_ac_upgene$V1, vpa_ac_upgene$V1)
co_ac_downgene <- intersect(tsa_ac_downgene$V1, vpa_ac_downgene$V1)
co_vg_upgene <- intersect(tsa_vg_upgene$V1, vpa_vg_upgene$V1)
co_vg_downgene <- intersect(tsa_vg_downgene$V1, vpa_vg_downgene$V1)

#####Look at all these genes as a whole######
co_de_genes <- c(co_ac_upgene, co_vg_upgene, co_ac_downgene, co_vg_downgene)
co_de_genes <- unique(co_de_genes)

####import Hdac1 classes of genes####
class1_gene <- read.delim("~/GDrive files/foxh1 and hdac1/compare VPA and TSA/Hdac1 class DE genes/class1_gene_both_h3k27.txt", header = F)
class2_gene <- read.delim("~/GDrive files/foxh1 and hdac1/compare VPA and TSA/Hdac1 class DE genes/class2_gene_h3k27me3_only.txt", header = F)
class3_gene <- read.delim("~/GDrive files/foxh1 and hdac1/compare VPA and TSA/Hdac1 class DE genes/class3_gene_h3k27ac_only.txt", header = F)

class1_degene <- subset(class1_gene, class1_gene$V1 %in% co_de_genes)
class2_degene <- subset(class2_gene, class2_gene$V1 %in% co_de_genes)
class3_degene <- subset(class3_gene, class3_gene$V1 %in% co_de_genes)

####import TPMs for heatmaps####
tsa_ac_tpm <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/ectoderm_TPM_DMSO_vs_TSA.txt", header=T)
tsa_ac_tpm[, 1] <- str_remove(tsa_ac_tpm$gene_id, "gene-")
tsa_vg_tpm <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/endoderm_TPM_DMSO_vs_TSA.txt", header=T)
tsa_vg_tpm[, 1] <- str_remove(tsa_vg_tpm$gene_id, "gene-")
vpa_ac_tpm <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/ectoderm_TPM_h2o_vs_VPA.txt", header=T)
vpa_ac_tpm[, 1] <- str_remove(vpa_ac_tpm$gene_id, "gene-")
vpa_vg_tpm <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/endoderm_TPM_h2o_vs_VPA.txt", header=T)
vpa_vg_tpm[, 1] <- str_remove(vpa_vg_tpm$gene_id, "gene-")

####Sum up TPM in AC and VG for heatmap plotting as a whole####
tsa_tpm <- data.frame(matrix(nrow= nrow(tsa_ac_tpm), ncol= 5))
tsa_tpm[, 1] <- tsa_ac_tpm$gene_id
tsa_tpm[, 2] <- tsa_ac_tpm$TPM + tsa_vg_tpm$TPM
tsa_tpm[, 3] <- tsa_ac_tpm$TPM.1 + tsa_vg_tpm$TPM.1
tsa_tpm[, 4] <- tsa_ac_tpm$TPM.2 + tsa_vg_tpm$TPM.2
tsa_tpm[, 5] <- tsa_ac_tpm$TPM.3 + tsa_vg_tpm$TPM.3
colnames(tsa_tpm) <- c("gene", "dmso1", "dmso2", "tsa1", "tsa2")

vpa_tpm <- data.frame(matrix(nrow= nrow(vpa_ac_tpm), ncol= 5))
vpa_tpm[, 1] <- vpa_ac_tpm$gene_id
vpa_tpm[, 2] <- vpa_ac_tpm$TPM + vpa_vg_tpm$TPM
vpa_tpm[, 3] <- vpa_ac_tpm$TPM.1 + vpa_vg_tpm$TPM.1
vpa_tpm[, 4] <- vpa_ac_tpm$TPM.2 + vpa_vg_tpm$TPM.2
vpa_tpm[, 5] <- vpa_ac_tpm$TPM.3 + vpa_vg_tpm$TPM.3
colnames(vpa_tpm) <- c("gene", "h2o1", "h2o2", "vpa1", "vpa2")

tsa_vpa_tpm <- cbind(tsa_tpm, vpa_tpm[, 2:5])
###############################################################
####Class1 DEgene heatmap####
class1_degene_tpm <- subset(tsa_vpa_tpm, tsa_vpa_tpm$gene %in% class1_degene$V1)
class1_hp <- pheatmap(log(class1_degene_tpm[, 2:9]+1, 2), 
                     scale = "row",
                     cutree_rows = 2, 
                     cluster_cols = F, 
                     show_rownames = F,
                     clustering_method = "ward.D2", 
                     border_color = NA,
                     cellwidth = 20,
                     cellheight = 1,
                     treeheight_row = 0,
                     gaps_col = 4,
                     col=colorRampPalette(rev(brewer.pal(11,"BrBG")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/Hdac1 class DE genes/hp_Class1_deGenes.tiff', class1_hp, device = "tiff", dpi = 300)

class1_cluster <- cutree(class1_hp$tree_row, k=2)
sum(class1_cluster == 1)
sum(class1_cluster == 2)

####Class2 DEgene heatmap####
class2_degene_tpm <- subset(tsa_vpa_tpm, tsa_vpa_tpm$gene %in% class2_degene$V1)
class2_hp <- pheatmap(log(class2_degene_tpm[, 2:9]+1, 2), 
                      scale = "row",
                      cutree_rows = 2, 
                      cluster_cols = F, 
                      show_rownames = F,
                      clustering_method = "ward.D2", 
                      border_color = NA,
                      cellwidth = 20,
                      cellheight = 1,
                      treeheight_row = 0,
                      gaps_col = 4,
                      col=colorRampPalette(rev(brewer.pal(11,"BrBG")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/Hdac1 class DE genes/hp_Class2_deGenes.tiff', class2_hp, device = "tiff", dpi = 300)

class2_cluster <- cutree(class2_hp$tree_row, k=2)
sum(class2_cluster == 1)
sum(class2_cluster == 2)

####Class3 DEgene heatmap####
class3_degene_tpm <- subset(tsa_vpa_tpm, tsa_vpa_tpm$gene %in% class3_degene$V1)
class3_hp <- pheatmap(log(class3_degene_tpm[, 2:9]+1, 2), 
                      scale = "row",
                      cutree_rows = 2, 
                      cluster_cols = F, 
                      show_rownames = F,
                      clustering_method = "ward.D2", 
                      border_color = NA,
                      cellwidth = 20,
                      cellheight = 1,
                      treeheight_row = 0,
                      gaps_col = 4,
                      col=colorRampPalette(rev(brewer.pal(11,"BrBG")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/Hdac1 class DE genes/hp_Class3_deGenes.tiff', class3_hp, device = "tiff", dpi = 300)

class3_cluster <- cutree(class3_hp$tree_row, k=2)
sum(class3_cluster == 1)
sum(class3_cluster == 2)















