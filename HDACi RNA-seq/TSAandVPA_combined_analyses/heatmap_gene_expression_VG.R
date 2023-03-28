################################################################################
##################Compare TSA DEgenes and VPA DEgenes in VG#####################
################################################################################
tsa_degene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/VG/DEgene_TSA_VG.txt", header=T)
vpa_degene <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/VG/DEgene_VPA_VG.txt", header=T)

tsa_tpmTable <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/VG/endoderm_TPM_DMSO_vs_TSA.txt", header=T)
vpa_tpmTable <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/VG/endoderm_TPM_h2o_vs_VPA.txt", header=T)

#######Set up 3 gene groups######
tsa_only_gene <- subset(tsa_degene, !(rownames(tsa_degene) %in% rownames(vpa_degene)))
tsa_only_gene <- rownames(tsa_only_gene)

shared_gene <- subset(tsa_degene, rownames(tsa_degene) %in% rownames(vpa_degene))
shared_gene <- rownames(shared_gene)

vpa_only_gene <- subset(vpa_degene, !(rownames(vpa_degene) %in% rownames(tsa_degene)))
vpa_only_gene <- rownames(vpa_only_gene)

########Extract TPM values for these genes######
tsa_only_gene_tpm1 <- subset(tsa_tpmTable, tsa_tpmTable$gene_id %in% tsa_only_gene)
tsa_only_gene_tpm2 <- subset(vpa_tpmTable, vpa_tpmTable$gene_id %in% tsa_only_gene)
tsa_only_gene_tpm <- cbind(tsa_only_gene_tpm1, tsa_only_gene_tpm2)
col <- c(1:5, 7:10)
tsa_only_gene_tpm <- tsa_only_gene_tpm[, col]
colnames(tsa_only_gene_tpm) <- c("geneID", "dmso1", "dmso2", "tsa1", "tsa2",
                                 "h2o1", "h2o2", "vpa1", "vpa2")
rownames(tsa_only_gene_tpm) <- tsa_only_gene_tpm[, 1]
tsa_only_gene_tpm <- tsa_only_gene_tpm[, -1]

shared_gene_tpm1 <- subset(tsa_tpmTable, tsa_tpmTable$gene_id %in% shared_gene)
shared_gene_tpm2 <- subset(vpa_tpmTable, vpa_tpmTable$gene_id %in% shared_gene)
shared_gene_tpm <- cbind(shared_gene_tpm1, shared_gene_tpm2)
col <- c(1:5, 7:10)
shared_gene_tpm <- shared_gene_tpm[, col]
colnames(shared_gene_tpm) <- c("geneID", "dmso1", "dmso2", "tsa1", "tsa2",
                               "h2o1", "h2o2", "vpa1", "vpa2")
rownames(shared_gene_tpm) <- shared_gene_tpm[, 1]
shared_gene_tpm <- shared_gene_tpm[, -1]

vpa_only_gene_tpm1 <- subset(tsa_tpmTable, tsa_tpmTable$gene_id %in% vpa_only_gene)
vpa_only_gene_tpm2 <- subset(vpa_tpmTable, vpa_tpmTable$gene_id %in% vpa_only_gene)
vpa_only_gene_tpm <- cbind(vpa_only_gene_tpm1, vpa_only_gene_tpm2)
col <- c(1:5, 7:10)
vpa_only_gene_tpm <- vpa_only_gene_tpm[, col]
colnames(vpa_only_gene_tpm) <- c("geneID", "dmso1", "dmso2", "tsa1", "tsa2",
                                 "h2o1", "h2o2", "vpa1", "vpa2")
rownames(vpa_only_gene_tpm) <- vpa_only_gene_tpm[, 1]
vpa_only_gene_tpm <- vpa_only_gene_tpm[, -1]


############Heatmaps#############
library(pheatmap)
library(ggplot2)
library("RColorBrewer")
tsa_only_hp <- pheatmap(log(tsa_only_gene_tpm[, 1:4]+1, 2), 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        clustering_method = "ward.D2", 
                        border_color = NA,
                        cellwidth = 10,
                        cellheight = 0.2,
                        treeheight_row = 0,
                        col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255))
reordered_tsa_only_gene_tpm <- tsa_only_gene_tpm[tsa_only_hp$tree_row$order, ]
tsa_only_hp <- pheatmap(log(reordered_tsa_only_gene_tpm+1, 2), 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        cluster_rows = F, 
                        border_color = NA,
                        cellwidth = 10,
                        cellheight = 0.2,
                        treeheight_row = 0,
                        gaps_col = 4,
                        col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/VG/hp_tsa_only_Degenes.tiff', tsa_only_hp, device = "tiff", dpi = 300)

################################################################################
shared_hp <- pheatmap(log(shared_gene_tpm+1, 2), 
                      scale = "row",
                      cutree_rows = 1, 
                      cluster_cols = F, 
                      show_rownames = F,
                      clustering_method = "ward.D2", 
                      border_color = NA,
                      cellwidth = 10,
                      cellheight = 0.2,
                      treeheight_row = 0,
                      gaps_col = 4,
                      col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/VG/hp_shared_Degenes.tiff', shared_hp, device = "tiff", dpi = 300)

################################################################################
vpa_only_hp <- pheatmap(log(vpa_only_gene_tpm[, 5:8]+1, 2), 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        clustering_method = "ward.D2", 
                        border_color = NA,
                        cellwidth = 10,
                        cellheight = 0.2,
                        treeheight_row = 0,
                        gaps_col = 4,
                        col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255))
reordered_vpa_only_gene_tpm <- vpa_only_gene_tpm[vpa_only_hp$tree_row$order, ]
vpa_only_hp <- pheatmap(log(reordered_vpa_only_gene_tpm+1, 2), 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        cluster_rows = F, 
                        border_color = NA,
                        cellwidth = 10,
                        cellheight = 0.2,
                        treeheight_row = 0,
                        gaps_col = 4,
                        col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/VG/hp_vpa_only_Degenes.tiff', vpa_only_hp, device = "tiff", dpi = 300)





