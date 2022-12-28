##############st8################
st8_annotate_output <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/distribution/st8_hdac1_wt_IDR0.05_in_merged_peaks.annotate", header=FALSE)
st8 <- data.frame(st8_annotate_output$V8)
st8 <- data.frame(st8[-1,])
colnames(st8) <- c("region")
st8$region <- as.character(st8$region)
st8 <- data.frame(sapply(strsplit(st8$region, "[()]"), "[", 1))
colnames(st8) <- c("region")
st8$region <- as.character(st8$region)
type <- unique(st8)
library(stringr)
exon <- sum(str_count(st8, "exon"))
intergenic <- sum(str_count(st8, "Intergenic"))
intron <- sum(str_count(st8, "intron"))
TSS <- sum(str_count(st8, "TSS"))
TTS <- sum(str_count(st8, "TTS"))
data <- data.frame(
  type= type[1:5,],
  value= c(571,480,145,114,27)
)
library(ggplot2)
# Basic piechart
pie_st8 <- ggplot(data, aes(x="", y=value, fill=type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
pie_st8
ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/distribution/st8_pie_on_IDR_peaks.tiff', pie_st8, device = "tiff", dpi = 100)



##############st9################
st9_annotate_output <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/distribution/st9_hdac1_wt_IDR0.05_in_merged_peaks.annotate", header=FALSE)
st9 <- data.frame(st9_annotate_output$V8)
st9 <- data.frame(st9[-1,])
colnames(st9) <- c("region")
st9$region <- as.character(st9$region)
st9 <- data.frame(sapply(strsplit(st9$region, "[()]"), "[", 1))
colnames(st9) <- c("region")
st9$region <- as.character(st9$region)
type <- unique(st9)
library(stringr)
exon <- sum(str_count(st9, "exon"))
intergenic <- sum(str_count(st9, "Intergenic"))
intron <- sum(str_count(st9, "intron"))
TSS <- sum(str_count(st9, "TSS"))
TTS <- sum(str_count(st9, "TTS"))
data <- data.frame(
  type= type[1:5,],
  value= c(1919,2093,338,137,648)
)
library(ggplot2)
# Basic piechart
pie_st9 <- ggplot(data, aes(x="", y=value, fill=type)) +
           geom_bar(stat="identity", width=1, color="white") +
           coord_polar("y", start=0) +
           theme_void() # remove background, grid, numeric labels
pie_st9
ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/distribution/st9_pie_on_IDR_peaks.tiff', pie_st9, device = "tiff", dpi = 100)





#############st10.5#################
st105_annotate_output <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/distribution/st105_hdac1_wt_IDR0.05_in_merged_peaks.annotate", header=FALSE)
st105 <- data.frame(st105_annotate_output$V8)
st105 <- data.frame(st105[-1,])
colnames(st105) <- c("region")
st105$region <- as.character(st105$region)
st105 <- data.frame(sapply(strsplit(st105$region, "[()]"), "[", 1))
colnames(st105) <- c("region")
st105$region <- as.character(st105$region)
type <- unique(st105)
library(stringr)
exon <- sum(str_count(st105, "exon"))
intergenic <- sum(str_count(st105, "Intergenic"))
intron <- sum(str_count(st105, "intron"))
TSS <- sum(str_count(st105, "TSS"))
TTS <- sum(str_count(st105, "TTS"))
data <- data.frame(
  type= type[1:5,],
  value= c(6343,8056,5525,2400,354)
)
library(ggplot2)
# Basic piechart
pie_st105 <- ggplot(data, aes(x="", y=value, fill=type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels
pie_st105
ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/distribution/st105_pie_on_IDR_peaks.tiff', pie_st105, device = "tiff", dpi = 100)



