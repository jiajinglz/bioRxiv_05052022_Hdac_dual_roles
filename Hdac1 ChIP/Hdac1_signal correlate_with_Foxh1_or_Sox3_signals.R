####################################################################
###################Signals on HDAC1 IDR peaks st9###################
####################################################################
counts <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 chip seq/overlap foxh1 sox3 at st9/counts_st9_hdac1_foxh1_sox3_signals_on_hdac1_peaks.bed", header=FALSE)
peak_length <- as.data.frame(0.001*(counts$V3-counts$V2))
counts <- cbind(counts, peak_length)
colnames(counts) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8')
signal_table <- as.data.frame(matrix(nrow = nrow(counts), ncol = 4))
signal_table[,1] <- paste(counts$V1, counts$V2, counts$V3, sep = "\t")
signal_table[,2] <- log2(counts$V5/ counts$V8+1)
signal_table[,3] <- log2(counts$V6/ counts$V8+1)
signal_table[,4] <- log2(counts$V7/ counts$V8+1)
colnames(signal_table) <- c('peak', 'hdac1', 'foxh1', 'sox3')
library(ggplot2)
library(ggpmisc)
plot <-
  ggplot(signal_table, aes(x=foxh1, y=hdac1))+
  geom_point(size= 0.5, color="black", alpha =1/10)+
  geom_smooth(method= 'lm', color= "tomato", formula = y ~ x)+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  theme_minimal()+
  theme(aspect.ratio = 1)
plot


plot <-
  ggplot(signal_table, aes(x=foxh1, y=hdac1))+
  geom_point(size= 0.5, color="black", alpha =1/10)+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method= 'lm', color= "tomato", formula = y ~ x)+
  stat_smooth_func(geom= "text", method= "lm", hjust=0, parse=T)+
  theme_minimal()+
  theme(aspect.ratio = 1)
plot











