st9IDRgenes <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/st9_hdac1_IDR_peak_closest_gene_body.bed", header=F)
nCounts <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/normalized_counts_timeCourse.txt", header=F)
which(is.na(match(st9IDRgenes$V4, nCounts$V1)))###quick check to make sure no N/A match
st9IDRgeneCounts <- nCounts[match(st9IDRgenes$V4, nCounts$V1),]
st9NonIDRgeneCounts <- subset(nCounts, !(nCounts$V1 %in% st9IDRgeneCounts$V1))
library(dplyr)
st9IDRgeneCounts <- st9IDRgeneCounts %>% mutate(type= c("IDR") )
st9NonIDRgeneCounts <- st9NonIDRgeneCounts %>% mutate(type= c("nonIDR") )
st9Counts <- rbind(st9IDRgeneCounts, st9NonIDRgeneCounts)
#######replace 0 with 1 for log transformation###########
st9Counts[st9Counts == 0] <- 1
logSt9Counts <- log2(st9Counts[,2:9])
col <- c(1,10)
logSt9Counts <- cbind(st9Counts[,col], logSt9Counts)
colnames(logSt9Counts) <- c("geneName", "type", "four1", "four2", "five1", "five2",
                            "fihalf1", "fihalf2", "seven1", "seven2")
############plot######################
library(ggplot2)
# Basic violin plot
p9 <- ggplot(logSt9Counts, aes(x= type, y= five1, fill= type))+ 
  geom_violin(trim = F, scale = "count")+
  stat_summary(fun=mean, geom="point", size=1, color= "red")+
  scale_fill_manual(values=c("salmon1", "skyblue1"))+
  theme_minimal()+
  theme(aspect.ratio = 1)
p9
ggsave('C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/st9_IDR_peak_on_5hpf.tiff', p9, device = "tiff", dpi = 500)
#############t-test--unpaired############
x<- st9IDRgeneCounts$V4
y<- st9NonIDRgeneCounts$V4
t9 <- t.test(x, y, alternative = "two.sided", var.equal = FALSE)
t9
#############splitViolin##############



################st10.5####################
st105IDRgenes <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/st105_hdac1_IDR_peak_closest_gene_body.bed", header=F)
nCounts <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/normalized_counts_timeCourse.txt", header=F)
which(is.na(match(st105IDRgenes$V4, nCounts$V1)))###quick check to make sure no N/A match
st105IDRgeneCounts <- nCounts[match(st105IDRgenes$V4, nCounts$V1),]
st105NonIDRgeneCounts <- subset(nCounts, !(nCounts$V1 %in% st105IDRgeneCounts$V1))
library(dplyr)
st105IDRgeneCounts <- st105IDRgeneCounts %>% mutate(type= c("IDR") )
st105NonIDRgeneCounts <- st105NonIDRgeneCounts %>% mutate(type= c("nonIDR") )
st105Counts <- rbind(st105IDRgeneCounts, st105NonIDRgeneCounts)
#######replace 0 with 1 for log transformation###########
st105Counts[st105Counts == 0] <- 1
logSt105Counts <- log2(st105Counts[,2:9])
col <- c(1,10)
logSt105Counts <- cbind(st105Counts[,col], logSt105Counts)
colnames(logSt105Counts) <- c("geneName", "type", "four1", "four2", "five1", "five2",
                            "fihalf1", "fihalf2", "seven1", "seven2")
############plot######################
library(ggplot2)
# Basic violin plot
p105 <- ggplot(logSt105Counts, aes(x= type, y= seven1, fill= type))+ 
  geom_violin(trim = F, scale = "count")+
  stat_summary(fun=mean, geom="point", size=1, color= "red")+
  scale_fill_manual(values=c("salmon1", "skyblue1"))+
  theme_minimal()+
  theme(aspect.ratio = 1)
p105
ggsave('C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/binding correlate expression/st105_IDR_peak_on_7hpf.tiff', p105, device = "tiff", dpi = 500)
#############t-test--unpaired############
x<- st105IDRgeneCounts$V8
y<- st105NonIDRgeneCounts$V8
t105 <- t.test(x, y, alternative = "two.sided", var.equal = FALSE)
t105
t105$p.value




