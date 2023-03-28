#####################################################################
####################Heatmap for Category#############################
#####################################################################
ds_tpm <- read.delim("G:/My Drive/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/TPM_table_ira_dissected_st105_tissues.txt", header=FALSE)
convert <- read.delim("G:/My Drive/foxh1 and hdac1/XENTR_10.0_GeneID_GeneName_converter.txt", header=F)
geneName <- convert[match(ds_tpm$V1, convert$V2),4]
ds_tpm <- cbind(geneName, ds_tpm[, 2:11])
colnames(ds_tpm) <- c("geneName", "ac1", "ac2", "dm1", "dm2", "lm1", "lm2",
                      "vm1", "vm2", "vg1", "vg2")
ds_tpm <- ds_tpm[complete.cases(ds_tpm),]
##########make a new data frame from mean TPM since rep variation is too large######
ds_tpm1 <- data.frame(matrix(nrow = nrow(ds_tpm), ncol=4))
ds_tpm1[, 1] <- ds_tpm$geneName
ds_tpm1[, 2] <- rowMeans(ds_tpm[2:3])
ds_tpm1[, 3] <- rowMeans(ds_tpm[4:9])
ds_tpm1[, 4] <- rowMeans(ds_tpm[10:11])
colnames(ds_tpm1) <- c("geneName", "ac", "mz", "vg")
##########find lowly expressed genes by TPM< 1 in any given tissues#########
library(rje)
low_expression_genes <- subset(ds_tpm1, rowMaxs(ds_tpm1[, 2:4])< 1)
ds_tpm2 <- subset(ds_tpm1, rowMaxs(ds_tpm1[, 2:4])>= 1)
#write.table(low_expression_genes, file= 'G:/My Drive/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/low expression genes.txt', quote = F, row.names = F, sep = "\t")
##########find low variant genes by CV< 10% TPM ##############
low_variant_genes <- ds_tpm2[apply(ds_tpm2[, 2:4] , 1 , function(x) abs(sd(x)/mean(x))< 0.1 ), ]
ds_tpm3 <- ds_tpm2[apply(ds_tpm2[, 2:4] , 1 , function(x) abs(sd(x)/mean(x))>= 0.1 ), ]
#write.table(low_variant_genes, file= 'G:/My Drive/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/low variation genes.txt', quote = F, row.names = F, sep = "\t")
###########simple method that looks for max TPM in a tissue###########
rownames(ds_tpm3) <- ds_tpm3[, 1]
ds_tpm3 <- ds_tpm3[, -1]
tissue <- colnames(ds_tpm3)[apply(ds_tpm3, 1, which.max)]
ds_tpm4 <- cbind(ds_tpm3, tissue)
ac_genes <- subset(ds_tpm4, ds_tpm4$tissue == "ac")
#write.table(ac_genes, file= 'G:/My Drive/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/ac localized genes.txt', quote = F, row.names = T, sep = "\t")
mz_genes <- subset(ds_tpm4, ds_tpm4$tissue == "mz")
#write.table(mz_genes, file= 'G:/My Drive/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/mz localized genes.txt', quote = F, row.names = T, sep = "\t")
vg_genes <- subset(ds_tpm4, ds_tpm4$tissue == "vg")
#write.table(vg_genes, file= 'G:/My Drive/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/vg localized genes.txt', quote = F, row.names = T, sep = "\t")







##########TRASH BELOW############







###########hp for germ layer specific genes################
log_ds_tpm <- cbind(ds_tpm3$geneName, log(ds_tpm3[, 2:4]+1, 2))
log_ds_tpm <- log_ds_tpm[!apply(log_ds_tpm[, 2:4] , 1 , function(x) sd(x)==0 ), ]
####################heatmap#######################
library(pheatmap)
library("RColorBrewer")
hp <- pheatmap(         log_ds_tpm[, 2:4], 
                        scale = "row",
                        cutree_rows = 3, 
                        cluster_cols = F, 
                        show_rownames = F,
                        clustering_method = "mcquitty", 
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

ggsave('G:/My Drive/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/hp_all_genes.tiff', hp, device = "tiff", dpi = 500)




