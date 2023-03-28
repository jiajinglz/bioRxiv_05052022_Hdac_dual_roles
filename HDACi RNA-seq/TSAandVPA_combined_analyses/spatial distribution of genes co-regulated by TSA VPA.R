################################################################################
##############Spatial distribution of DEgenes by TSA and VPA####################
################################################################################
#LOAD PACKAGES#
library(stringr)
library("ggVennDiagram")
library(ggplot2)
library(ggalluvial)
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
write.table(co_ac_upgene,file="~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/vpa_tsa_AC_coUpDEgenes.txt", quote = FALSE, row.names = F, col.names = F)
co_ac_downgene <- intersect(tsa_ac_downgene$V1, vpa_ac_downgene$V1)
write.table(co_ac_downgene,file="~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/vpa_tsa_AC_coDownDEgenes.txt", quote = FALSE, row.names = F, col.names = F)
co_vg_upgene <- intersect(tsa_vg_upgene$V1, vpa_vg_upgene$V1)
write.table(co_vg_upgene,file="~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/vpa_tsa_VG_coUpDEgenes.txt", quote = FALSE, row.names = F, col.names = F)
co_vg_downgene <- intersect(tsa_vg_downgene$V1, vpa_vg_downgene$V1)
write.table(co_vg_downgene,file="~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/vpa_tsa_VG_coDownDEgenes.txt", quote = FALSE, row.names = F, col.names = F)

####Four-way pie####
mydata <- list(
  A= co_ac_upgene, 
  B= co_ac_downgene, 
  C= co_vg_upgene, 
  D= co_vg_downgene)

p <- ggVennDiagram(mydata, label_alpha = 0, label = "count")+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values=c("#1E90FF", "#1E90FF", "#DAA520", "#DAA520"))
p
ggsave('~/GDrive files/foxh1 and hdac1/compare VPA and TSA/Venn_DEgenes_in_TSA_VPA.tiff', p, device = "tiff", dpi = 300)

p <- ggVennDiagram(mydata, label_alpha = 0, label = "none")+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values=c("#F08080", "#F08080", "#3CB371", "#3CB371"))
p

ggsave('~/GDrive files/foxh1 and hdac1/compare VPA and TSA/Venn_DEgenes_in_TSA_VPA_no#.tiff', p, device = "tiff", dpi = 300)

##########################################
####Heatmap for selected example genes####
##########################################
spatial_expression <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/TPM_table_ira_dissected_st105_tissues.txt", header=F)
spatial_expression[, 1] <- str_remove(spatial_expression$V1, "gene-")
colnames(spatial_expression) <- c("gene", "ac1", "ac2", "dm1", "dm2", "lm1", "lm2",
                                  "vm1", "vm2", "vg1", "vg2")
co_ac_upgene_spatial <- subset(spatial_expression, spatial_expression$gene %in% co_ac_upgene)
co_ac_downgene_spatial <- subset(spatial_expression, spatial_expression$gene %in% co_ac_downgene)
co_vg_upgene_spatial <- subset(spatial_expression, spatial_expression$gene %in% co_vg_upgene)
co_vg_downgene_spatial <- subset(spatial_expression, spatial_expression$gene %in% co_vg_downgene)

ac_up_gene_eg <- c("adrb2", "foxc1", "hhex", "osr1", "snai1", "vegfa")
ac_down_gene_eg <- c("hes4", "bambi", "dlx5", "foxi4.2", "irx2", "tfap2a")
vg_up_gene_eg <- c("jun", "numbl", "esr-5", "cdx2")
vg_down_gene_eg <- c("foxa2", "gata5", "osr1", "pnhd")

tsa_ac_tpm <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/ectoderm_TPM_DMSO_vs_TSA.txt", header=T)
tsa_ac_tpm[, 1] <- str_remove(tsa_ac_tpm$gene_id, "gene-")
tsa_vg_tpm <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/endoderm_TPM_DMSO_vs_TSA.txt", header=T)
tsa_vg_tpm[, 1] <- str_remove(tsa_vg_tpm$gene_id, "gene-")
vpa_ac_tpm <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/ectoderm_TPM_h2o_vs_VPA.txt", header=T)
vpa_ac_tpm[, 1] <- str_remove(vpa_ac_tpm$gene_id, "gene-")
vpa_vg_tpm <- read.delim("~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/endoderm_TPM_h2o_vs_VPA.txt", header=T)
vpa_vg_tpm[, 1] <- str_remove(vpa_vg_tpm$gene_id, "gene-")

ac_up_gene_tsa_tpm <- subset(tsa_ac_tpm, tsa_ac_tpm$gene_id %in% ac_up_gene_eg)
ac_up_gene_vpa_tpm <- subset(vpa_ac_tpm, vpa_ac_tpm$gene_id %in% ac_up_gene_eg)
ac_up_gene_eg_tpm <- cbind(ac_up_gene_tsa_tpm, ac_up_gene_vpa_tpm[, 2:5])
colnames(ac_up_gene_eg_tpm) <- c("geneID", "dmso1", "dmso2", "tsa1", "tsa2",
                                 "h2o1", "h2o2", "vpa1", "vpa2")
rownames(ac_up_gene_eg_tpm) <- ac_up_gene_eg_tpm[, 1]
ac_up_gene_eg_tpm <- ac_up_gene_eg_tpm[, -1]
ac_up_hp <- pheatmap(log(ac_up_gene_eg_tpm+1, 2), 
                     scale = "row",
                     cutree_rows = 1, 
                     cluster_cols = F, 
                     show_rownames = T,
                     clustering_method = "ward.D2", 
                     border_color = "white",
                     cellwidth = 20,
                     cellheight = 20,
                     treeheight_row = 0,
                     gaps_col = 4,
                     col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/hp_ac_up_exampleGenes.tiff', ac_up_hp, device = "tiff", dpi = 300)


ac_down_gene_tsa_tpm <- subset(tsa_ac_tpm, tsa_ac_tpm$gene_id %in% ac_down_gene_eg)
ac_down_gene_vpa_tpm <- subset(vpa_ac_tpm, vpa_ac_tpm$gene_id %in% ac_down_gene_eg)
ac_down_gene_eg_tpm <- cbind(ac_down_gene_tsa_tpm, ac_down_gene_vpa_tpm[, 2:5])
colnames(ac_down_gene_eg_tpm) <- c("geneID", "dmso1", "dmso2", "tsa1", "tsa2",
                                   "h2o1", "h2o2", "vpa1", "vpa2")
rownames(ac_down_gene_eg_tpm) <- ac_down_gene_eg_tpm[, 1]
ac_down_gene_eg_tpm <- ac_down_gene_eg_tpm[, -1]
ac_down_hp <- pheatmap(log(ac_down_gene_eg_tpm+1, 2), 
                       scale = "row",
                       cutree_rows = 1, 
                       cluster_cols = F, 
                       show_rownames = T,
                       clustering_method = "ward.D2", 
                       border_color = "white",
                       cellwidth = 20,
                       cellheight = 20,
                       treeheight_row = 0,
                       gaps_col = 4,
                       col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/hp_ac_down_exampleGenes.tiff', ac_down_hp, device = "tiff", dpi = 300)


vg_up_gene_tsa_tpm <- subset(tsa_vg_tpm, tsa_ac_tpm$gene_id %in% vg_up_gene_eg)
vg_up_gene_vpa_tpm <- subset(vpa_vg_tpm, vpa_ac_tpm$gene_id %in% vg_up_gene_eg)
vg_up_gene_eg_tpm <- cbind(vg_up_gene_tsa_tpm, vg_up_gene_vpa_tpm[, 2:5])
colnames(vg_up_gene_eg_tpm) <- c("geneID", "dmso1", "dmso2", "tsa1", "tsa2",
                                 "h2o1", "h2o2", "vpa1", "vpa2")
rownames(vg_up_gene_eg_tpm) <- vg_up_gene_eg_tpm[, 1]
vg_up_gene_eg_tpm <- vg_up_gene_eg_tpm[, -1]
vg_up_hp <- pheatmap(log(vg_up_gene_eg_tpm+1, 2), 
                     scale = "row",
                     cutree_rows = 1, 
                     cluster_cols = F, 
                     show_rownames = T,
                     clustering_method = "ward.D2", 
                     border_color = "white",
                     cellwidth = 20,
                     cellheight = 20,
                     treeheight_row = 0,
                     gaps_col = 4,
                     col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/hp_vg_up_exampleGenes.tiff', vg_up_hp, device = "tiff", dpi = 300)


vg_down_gene_tsa_tpm <- subset(tsa_vg_tpm, tsa_ac_tpm$gene_id %in% vg_down_gene_eg)
vg_down_gene_vpa_tpm <- subset(vpa_vg_tpm, vpa_ac_tpm$gene_id %in% vg_down_gene_eg)
vg_down_gene_eg_tpm <- cbind(vg_down_gene_tsa_tpm, vg_down_gene_vpa_tpm[, 2:5])
colnames(vg_down_gene_eg_tpm) <- c("geneID", "dmso1", "dmso2", "tsa1", "tsa2",
                                   "h2o1", "h2o2", "vpa1", "vpa2")
rownames(vg_down_gene_eg_tpm) <- vg_down_gene_eg_tpm[, 1]
vg_down_gene_eg_tpm <- vg_down_gene_eg_tpm[, -1]
vg_down_hp <- pheatmap(log(vg_down_gene_eg_tpm+1, 2), 
                       scale = "row",
                       cutree_rows = 1, 
                       cluster_cols = F, 
                       show_rownames = T,
                       clustering_method = "ward.D2", 
                       border_color = "white",
                       cellwidth = 20,
                       cellheight = 20,
                       treeheight_row = 0,
                       gaps_col = 4,
                       col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(255))
ggsave('~/GDrive files/foxh1 and hdac1/Compare VPA and TSA/germ layer specific/hp_vg_down_exampleGenes.tiff', vg_down_hp, device = "tiff", dpi = 300)

################################################################################
################################################################################
################################################################################
###########Spatial RNA profile#############

####################################
#####Import Class 1 Hdac1 genes#####
####################################
class1_gene <- read.delim("~/GDrive files/foxh1 and hdac1/compare VPA and TSA/germ layer specific/class1_gene_both_h3k27.txt", header = F)

####Find DEgenes in Class 1 Hdac1 genes####
class1_ac_up <- intersect(class1_gene$V1, co_ac_upgene)
class1_ac_down <- intersect(class1_gene$V1, co_ac_downgene)
class1_vg_up <- intersect(class1_gene$V1, co_vg_upgene)
class1_vg_down <- intersect(class1_gene$V1, co_vg_downgene)

####Germ layer categories#####
no_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/low expression genes.txt', header = T)
no_expression <- data.frame(V1= no_expression$geneName)
even_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/low variation genes.txt', header = T)
even_expression <- data.frame(V1= even_expression$geneName)
ac_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/ac localized genes.txt', header = T)
ac_expression <- data.frame(V1= rownames(ac_expression))
mz_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/mz localized genes.txt', header = T)
mz_expression <- data.frame(V1= rownames(mz_expression))
vg_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/vg localized genes.txt', header = T)
vg_expression <- data.frame(V1= rownames(vg_expression))

####assign Class1 co-DEgenes to spatial catergories####
class1_ac_up_ecto <- intersect(class1_ac_up, ac_expression$V1)
class1_ac_up_meso <- intersect(class1_ac_up, mz_expression$V1)
class1_ac_up_endo <- intersect(class1_ac_up, vg_expression$V1)
class1_ac_up_even <- intersect(class1_ac_up, even_expression$V1)
class1_ac_up_zero <- intersect(class1_ac_up, no_expression$V1)

class1_ac_down_ecto <- intersect(class1_ac_down, ac_expression$V1)
class1_ac_down_meso <- intersect(class1_ac_down, mz_expression$V1)
class1_ac_down_endo <- intersect(class1_ac_down, vg_expression$V1)
class1_ac_down_even <- intersect(class1_ac_down, even_expression$V1)
class1_ac_down_zero <- intersect(class1_ac_down, no_expression$V1)

class1_vg_up_ecto <- intersect(class1_vg_up, ac_expression$V1)
class1_vg_up_meso <- intersect(class1_vg_up, mz_expression$V1)
class1_vg_up_endo <- intersect(class1_vg_up, vg_expression$V1)
class1_vg_up_even <- intersect(class1_vg_up, even_expression$V1)
class1_vg_up_zero <- intersect(class1_vg_up, no_expression$V1)

class1_vg_down_ecto <- intersect(class1_vg_down, ac_expression$V1)
class1_vg_down_meso <- intersect(class1_vg_down, mz_expression$V1)
class1_vg_down_endo <- intersect(class1_vg_down, vg_expression$V1)
class1_vg_down_even <- intersect(class1_vg_down, even_expression$V1)
class1_vg_down_zero <- intersect(class1_vg_down, no_expression$V1)

##########summary plot table############
plot_table <- data.frame(matrix(nrow = 25, ncol=3))
plot_table[, 1] <- c("all", "all", "all", "all", "all", 
                     "ACup", "ACup", "ACup", "ACup", "ACup",
                     "ACdown", "ACdown", "ACdown", "ACdown", "ACdown",
                     "VGup", "VGup", "VGup", "VGup", "VGup",
                     "VGdown", "VGdown", "VGdown", "VGdown", "VGdown")
plot_table[, 2] <- c("AC", "MZ", "VG", "NL", "NP",
                     "AC", "MZ", "VG", "NL", "NP",
                     "AC", "MZ", "VG", "NL", "NP",
                     "AC", "MZ", "VG", "NL", "NP",
                     "AC", "MZ", "VG", "NL", "NP")
plot_table[, 3] <- c(length(ac_expression$V1), length(mz_expression$V1), length(vg_expression$V1),
                     length(even_expression$V1), length(no_expression$V1), length(class1_ac_up_ecto),
                     length(class1_ac_up_meso), length(class1_ac_up_endo), length(class1_ac_up_even),
                     length(class1_ac_up_zero), length(class1_ac_down_ecto), length(class1_ac_down_meso), 
                     length(class1_ac_down_endo), length(class1_ac_down_even), length(class1_ac_down_zero), 
                     length(class1_vg_up_ecto), length(class1_vg_up_meso), length(class1_vg_up_endo), 
                     length(class1_vg_up_even), length(class1_vg_up_zero), length(class1_vg_down_ecto), 
                     length(class1_vg_down_meso), length(class1_vg_down_endo), length(class1_vg_down_even), 
                     length(class1_vg_down_zero))
colnames(plot_table) <- c("group", "region", "frequency")


########Plot: Class 1 DEgenes##########
plot_table$region <- factor(plot_table$region, levels= c("AC", "MZ", "VG", "NL", "NP"))
plot_table$group <- factor(plot_table$group, levels= c("all", "ACup", "ACdown", "VGup", "VGdown"))
p <- ggplot(plot_table, aes(x= group, y= frequency, fill= region, width= 0.6))+
  geom_bar(position= "fill", stat= "identity")+
  scale_fill_manual(values=c("#F08080", "#F0E68C", "#3CB371", "#ADD8E6", "#C0C0C0"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 1/2)
p

ggsave('~/GDrive files/foxh1 and hdac1/compare VPA and TSA/germ layer specific/Class1_co_TSA_VPA_DEgenes_spatial_proportions.tiff', p, device = "tiff", dpi = 300)

#######Fisher's exact test for proportion significance#######
###number is success vs non-success in each col###
###AC###
fishertable <- matrix(c(9, 92-9, 3036, 28850-3036), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(34, 63-34, 3036, 28850-3036), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(8, 27-8, 3036, 28850-3036), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(10, 56-10, 3036, 28850-3036), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

###MZ###
fishertable <- matrix(c(7, 92-7, 965, 28850-965), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(14, 63-14, 965, 28850-965), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(3, 27-3, 965, 28850-965), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(10, 56-10, 965, 28850-965), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

###VG###
fishertable <- matrix(c(43, 92-43, 6170, 28850-6170), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(12, 63-12, 6170, 28850-6170), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(4, 27-4, 6170, 28850-6170), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(31, 56-31, 6170, 28850-6170), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft$p.value






####################################
#####Import Class 3 Hdac1 genes#####
####################################
class3_gene <- read.delim("~/GDrive files/foxh1 and hdac1/compare VPA and TSA/germ layer specific/class3_gene_h3k27ac_only.txt", header = F)

####Find DEgenes in Class 3 Hdac1 genes####
class3_ac_up <- intersect(class3_gene$V1, co_ac_upgene)
class3_ac_down <- intersect(class3_gene$V1, co_ac_downgene)
class3_vg_up <- intersect(class3_gene$V1, co_vg_upgene)
class3_vg_down <- intersect(class3_gene$V1, co_vg_downgene)

####Germ layer categories#####
no_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/low expression genes.txt', header = T)
no_expression <- data.frame(V1= no_expression$geneName)
even_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/low variation genes.txt', header = T)
even_expression <- data.frame(V1= even_expression$geneName)
ac_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/ac localized genes.txt', header = T)
ac_expression <- data.frame(V1= rownames(ac_expression))
mz_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/mz localized genes.txt', header = T)
mz_expression <- data.frame(V1= rownames(mz_expression))
vg_expression <- read.delim('~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ira dissection data heatmap to bar/vg localized genes.txt', header = T)
vg_expression <- data.frame(V1= rownames(vg_expression))

####assign Class3 co-DEgenes to spatial catergories####
class3_ac_up_ecto <- intersect(class3_ac_up, ac_expression$V1)
class3_ac_up_meso <- intersect(class3_ac_up, mz_expression$V1)
class3_ac_up_endo <- intersect(class3_ac_up, vg_expression$V1)
class3_ac_up_even <- intersect(class3_ac_up, even_expression$V1)
class3_ac_up_zero <- intersect(class3_ac_up, no_expression$V1)

class3_ac_down_ecto <- intersect(class3_ac_down, ac_expression$V1)
class3_ac_down_meso <- intersect(class3_ac_down, mz_expression$V1)
class3_ac_down_endo <- intersect(class3_ac_down, vg_expression$V1)
class3_ac_down_even <- intersect(class3_ac_down, even_expression$V1)
class3_ac_down_zero <- intersect(class3_ac_down, no_expression$V1)

class3_vg_up_ecto <- intersect(class3_vg_up, ac_expression$V1)
class3_vg_up_meso <- intersect(class3_vg_up, mz_expression$V1)
class3_vg_up_endo <- intersect(class3_vg_up, vg_expression$V1)
class3_vg_up_even <- intersect(class3_vg_up, even_expression$V1)
class3_vg_up_zero <- intersect(class3_vg_up, no_expression$V1)

class3_vg_down_ecto <- intersect(class3_vg_down, ac_expression$V1)
class3_vg_down_meso <- intersect(class3_vg_down, mz_expression$V1)
class3_vg_down_endo <- intersect(class3_vg_down, vg_expression$V1)
class3_vg_down_even <- intersect(class3_vg_down, even_expression$V1)
class3_vg_down_zero <- intersect(class3_vg_down, no_expression$V1)

##########summary plot table############
plot_table <- data.frame(matrix(nrow = 25, ncol=3))
plot_table[, 1] <- c("all", "all", "all", "all", "all", 
                     "ACup", "ACup", "ACup", "ACup", "ACup",
                     "ACdown", "ACdown", "ACdown", "ACdown", "ACdown",
                     "VGup", "VGup", "VGup", "VGup", "VGup",
                     "VGdown", "VGdown", "VGdown", "VGdown", "VGdown")
plot_table[, 2] <- c("AC", "MZ", "VG", "NL", "NP",
                     "AC", "MZ", "VG", "NL", "NP",
                     "AC", "MZ", "VG", "NL", "NP",
                     "AC", "MZ", "VG", "NL", "NP",
                     "AC", "MZ", "VG", "NL", "NP")
plot_table[, 3] <- c(length(ac_expression$V1), length(mz_expression$V1), length(vg_expression$V1),
                     length(even_expression$V1), length(no_expression$V1), length(class3_ac_up_ecto),
                     length(class3_ac_up_meso), length(class3_ac_up_endo), length(class3_ac_up_even),
                     length(class3_ac_up_zero), length(class3_ac_down_ecto), length(class3_ac_down_meso), 
                     length(class3_ac_down_endo), length(class3_ac_down_even), length(class3_ac_down_zero), 
                     length(class3_vg_up_ecto), length(class3_vg_up_meso), length(class3_vg_up_endo), 
                     length(class3_vg_up_even), length(class3_vg_up_zero), length(class3_vg_down_ecto), 
                     length(class3_vg_down_meso), length(class3_vg_down_endo), length(class3_vg_down_even), 
                     length(class3_vg_down_zero))
colnames(plot_table) <- c("group", "region", "frequency")


########Plot: Class 3 DEgenes##########
plot_table$region <- factor(plot_table$region, levels= c("AC", "MZ", "VG", "NL", "NP"))
plot_table$group <- factor(plot_table$group, levels= c("all", "ACup", "ACdown", "VGup", "VGdown"))
p <- ggplot(plot_table, aes(x= group, y= frequency, fill= region, width= 0.6))+
  geom_bar(position= "fill", stat= "identity")+
  scale_fill_manual(values=c("#F08080", "#F0E68C", "#3CB371", "#ADD8E6", "#C0C0C0"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 1/2)
p

ggsave('~/GDrive files/foxh1 and hdac1/compare VPA and TSA/germ layer specific/Class3_co_TSA_VPA_DEgenes_spatial_proportions.tiff', p, device = "tiff", dpi = 300)

#######Fisher's exact test for proportion significance#######
###number is success vs non-success in each col###
###AC###
fishertable <- matrix(c(17, 77-17, 3036, 28850-3036), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(51, 71-51, 3036, 28850-3036), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(2, 12-2, 3036, 28850-3036), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(16, 36-16, 3036, 28850-3036), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

###MZ###
fishertable <- matrix(c(10, 77-10, 965, 28850-965), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(7, 71-7, 965, 28850-965), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(3, 12-3, 965, 28850-965), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(8, 36-8, 965, 28850-965), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

###VG###
fishertable <- matrix(c(35, 77-35, 6170, 28850-6170), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(5, 71-5, 6170, 28850-6170), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(3, 12-3, 6170, 28850-6170), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft

fishertable <- matrix(c(9, 36-9, 6170, 28850-6170), nrow = 2)
ft <- fisher.test(fishertable, alternative = "greater")
ft






