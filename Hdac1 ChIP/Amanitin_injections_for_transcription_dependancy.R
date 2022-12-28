#################################################################################
############Compare Hdac1 enrichment between UT and a-Amanitin samples###########
#################################################################################
counts_table <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/amanitin/countsPerRegion_amanitinSamples.txt", header=T)
peak_length <- as.data.frame(0.001*(counts_table$X.end.-counts_table$X.start.))
colnames(peak_length) <- c('kb')
nrow(peak_length[peak_length$kb<0,]) ##check if any weird negative value
##wt= 56.785047 million reads
##ut= 58.984130 million reads
##a-amanitin= 65.668174 million reads
########find reads(counts) to peaks only#########
wtReads <- sum(counts_table$X.sampleWT.)
utReads <- sum(counts_table$X.sampleUT.)
amReads <- sum(counts_table$X.sampleAM.)
counts_table1 <- cbind(counts_table, peak_length)
##peak signal= counts/(peak kp* read depth), where read depth is # million reads mapped to peaks
##so peak signal is defined as counts per kilo base pairs per million reads mapped to peaks
wt_peak_signal <- as.data.frame(log2(counts_table1$X.sampleWT./(wtReads*10e-6*counts_table1$kb)+1))
ut_peak_signal <- as.data.frame(log2(counts_table1$X.sampleUT./(utReads*10e-6*counts_table1$kb)+1))
am_peak_signal <- as.data.frame(log2(counts_table1$X.sampleAM./(amReads*10e-6*counts_table1$kb)+1))
peak_signal_table <- cbind(counts_table[1:4], wt_peak_signal, ut_peak_signal, am_peak_signal)
colnames(peak_signal_table) <- c('chrm', 'start', 'end', 'name', 'wt', 'ut', 'amanitin')
#####quick t.test#####
wt_ut <- t.test(peak_signal_table$wt, peak_signal_table$ut)
wt_ut$p.value
ut_amanitin <- t.test(peak_signal_table$ut, peak_signal_table$amanitin)
ut_amanitin$p.value
wt_amanitin <- t.test(peak_signal_table$wt, peak_signal_table$amanitin)
wt_amanitin$p.value
###########box plot###############
library(dplyr)
wt_log_ncounts <- wt_peak_signal %>% mutate(type= c("WT") )
colnames(wt_log_ncounts) <- c('signal')
ut_log_ncounts <- ut_peak_signal %>% mutate(type= c("UT") )
colnames(ut_log_ncounts) <- c('signal')
am_log_ncounts <- am_peak_signal %>% mutate(type= c("a-Amanitin") )
colnames(am_log_ncounts) <- c('signal')
bp_table <- rbind(wt_log_ncounts, ut_log_ncounts, am_log_ncounts)
colnames(bp_table) <- c('signal', 'type')
bp_table$type <- factor(bp_table$type, levels=c("WT", "UT", "a-Amanitin"))
library(ggplot2)
# Basic box plot
bp <- ggplot(bp_table, aes(x= type, y= signal, fill= type))+ 
  geom_boxplot()+
  scale_fill_manual(values=c("azure3", "azure3", "azure3"))+
  theme_minimal()+
  theme(aspect.ratio = 1)+
  theme(axis.text = element_text(color="black"), text = element_text(size=9))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title = element_blank())
bp
ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/amanitin/boxplot_wt_ut_amanitin.tiff', bp, device = "tiff", dpi =200)
##################plot without WT#################
library(dplyr)
ut_log_ncounts <- ut_peak_signal %>% mutate(type= c("UT") )
colnames(ut_log_ncounts) <- c('signal')
am_log_ncounts <- am_peak_signal %>% mutate(type= c("a-Amanitin") )
colnames(am_log_ncounts) <- c('signal')
bp_table <- rbind(ut_log_ncounts, am_log_ncounts)
colnames(bp_table) <- c('signal', 'type')
bp_table$type <- factor(bp_table$type, levels=c("UT", "a-Amanitin"))
library(ggplot2)
# Basic box plot
bp <- ggplot(bp_table, aes(x= type, y= signal, fill= type))+ 
  geom_boxplot()+
  scale_fill_manual(values=c("azure3", "azure3"))+
  scale_y_continuous(limits = c(0,8))+
  theme_minimal()+
  theme(aspect.ratio = 1)+
  theme(axis.text = element_text(color="black"), text = element_text(size=12))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title = element_blank())
bp
ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 chip seq/amanitin/boxplot_ut_amanitin.tiff', bp, device = "tiff", dpi = 200)



















