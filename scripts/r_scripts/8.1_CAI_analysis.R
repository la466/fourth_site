# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/tiff/"), showWarnings = FALSE)

################
# File Inputs
################

compare_with_codonw <- read.csv("outputs/expression/_ecoli_codonw_compare.csv", head=TRUE)
cai_gc_file <- read.csv("outputs/expression/_CAI_a_ratios.csv", head=T)
CAI_a_nota_file <- read.csv("outputs/expression/_CAI_a_nota.csv", head = T)
CAI_high_low_file <- read.csv("outputs/expression/_high_low_CAI_compare.csv", head = T)





################
# Tests
################


sink("outputs/r_outputs/8.1_CAI_analysis.txt")


cai_max <- max(cai_gc_file$meanCAI)
cai_min <- min(cai_gc_file$meanCAI)

cat(sprintf('Max mean CAI value: %s\n', cai_max))
cat(sprintf('Min mean CAI value: %s\n', cai_min))


# Examine the correlation between the CAI values for each gene generated using
# own preferred codons and the deault values used by CodonW for Escherichia coli


cat("\n=================================================================================\n")
cat("Is there a strong correlation between CAI values calculated\n")
cat("using the 20 high expression reference set and default CodonW indices for E coli?\n")
cat("=================================================================================\n\n")

compare_cor_test <- cor.test(compare_with_codonw$my_cai, compare_with_codonw$codonw_cai, method = "spearman")
compare_cor_test



cat("\n=================================================================\n")
cat("Compare CAI values for genes with fourth site A and those without\n")
cat("=================================================================\n")

#########################################################################
# Comparing CAI values for those with fourth site A and those without
#########################################################################

# Compare the CAI values for genes with A and those without A

shapiro.test(CAI_a_nota_file$meanCAI_a)
shapiro.test(CAI_a_nota_file$meanCAI_not_a)


wilcox.test(CAI_a_nota_file$meanCAI_a, CAI_a_nota_file$meanCAI_not_a, paired = T)

mean_cai_a <- mean(CAI_a_nota_file$meanCAI_a)
mean_cai_nota <- mean(CAI_a_nota_file$meanCAI_not_a)
mean_diff <- mean(CAI_a_nota_file$meanCAI_a - CAI_a_nota_file$meanCAI_not_a)


cat(sprintf("Mean CAI +4A: %s +- %s\n", mean_cai_a, sd(CAI_a_nota_file$meanCAI_a)))
cat(sprintf("Mean CAI not +4A: %s +- %s\n", mean_cai_nota, sd(CAI_a_nota_file$meanCAI_not_a)))
cat(sprintf("Difference: %s\n\n", mean_diff))


# Test whether the mean CAIs in extreme GC genomes are significantly
# greater for genes with +4A

extremeGC3_a <- CAI_a_nota_file$meanCAI_a[CAI_a_nota_file$gc3>0.9 | CAI_a_nota_file$gc3<0.2]
extremeGC3_Nota <- CAI_a_nota_file$meanCAI_not_a[CAI_a_nota_file$gc3>0.9 | CAI_a_nota_file$gc3<0.2]
otherGC3_a <- CAI_a_nota_file$meanCAI_a[CAI_a_nota_file$gc3<=0.9 & CAI_a_nota_file$gc3>=0.2]
otherGC3_Nota <- CAI_a_nota_file$meanCAI_not_a[CAI_a_nota_file$gc3<=0.9 & CAI_a_nota_file$gc3>=0.2]







cat("\n=================================================================\n")
cat("Compare proportion +4 A in highly expressed genes and other genes\n")
cat("=================================================================\n")

shapiro.test(CAI_high_low_file$prop_a_high)
shapiro.test(CAI_high_low_file$prop_a_other)

wilcox.test(CAI_high_low_file$prop_a_high, CAI_high_low_file$prop_a_other, paired=T)

mean_aprop_high <- mean(CAI_high_low_file$prop_a_high)
mean_aprop_other <- mean(CAI_high_low_file$prop_a_other)
mean_prop_diff <- mean(CAI_high_low_file$prop_a_high - CAI_high_low_file$prop_a_other)



cat(sprintf("Mean proportion +4A high expressed genes: %s +- %s\n", mean_aprop_high, sd(CAI_high_low_file$prop_a_high)))
cat(sprintf("Mean proportion +4A other genes: %s +- %s\n", mean_aprop_other, sd(CAI_high_low_file$prop_a_other)))
cat(sprintf("Difference: %s\n\n", mean_prop_diff))



cat("\n=================================================================\n")
cat("mean CAI range of extreme GC3 genomes\n")
cat("=================================================================\n")

extreme_cai_max <- max(extremeGC3_a, extremeGC3_Nota)
extreme_cai_min <- min(extremeGC3_a, extremeGC3_Nota)
cat(sprintf('Max mean CAI value (gc3 extreme genomes): %s\n', extreme_cai_max))
cat(sprintf('Min mean CAI value (gc3 extreme genomes): %s\n', extreme_cai_min))

cat("\n===============================================\n")
cat("Compare between extreme GC genomes (< 0.2, > 0.9)\n")
cat("===============================================\n")

shapiro.test(extremeGC3_a)
shapiro.test(extremeGC3_Nota)
t.test(extremeGC3_a, extremeGC3_Nota, paired = T)

cat(sprintf("Mean CAI +4A extremes: %s +- %s\n", mean(extremeGC3_a), sd(extremeGC3_a)))
cat(sprintf("Mean CAI not +4A extremes: %s +- %s\n", mean(extremeGC3_Nota), sd(extremeGC3_Nota)))





highGCHighmeanCAI <- CAI_high_low_file$prop_a_high[CAI_high_low_file$gc3>0.9]
highGCLowmeanCAI <- CAI_high_low_file$prop_a_other[CAI_high_low_file$gc3>0.9]
lowGCHighmeanCAI <- CAI_high_low_file$prop_a_high[CAI_high_low_file$gc3<=0.2]
lowGCLowmeanCAI <- CAI_high_low_file$prop_a_other[CAI_high_low_file$gc3<=0.2]
extremesHigh <- CAI_high_low_file$prop_a_high[CAI_high_low_file$gc3>0.9 | CAI_high_low_file$gc3<0.2]
extremesOther <- CAI_high_low_file$prop_a_other[CAI_high_low_file$gc3>0.9 | CAI_high_low_file$gc3<0.2]
regularHigh <- CAI_high_low_file$prop_a_high[CAI_high_low_file$gc3<=0.9 & CAI_high_low_file$gc3>=0.2]
regularOther <- CAI_high_low_file$prop_a_other[CAI_high_low_file$gc3<=0.9 & CAI_high_low_file$gc3>=0.2]






cat("\n===============================================\n")
cat("Do extreme genomes have higher CAI for A genes?\n")
cat("===============================================\n")


lowGC <- CAI_a_nota_file[CAI_a_nota_file$gc3 < 0.2,]
highGC <- CAI_a_nota_file[CAI_a_nota_file$gc3 < 0.9,]

cat('Low GC\n')
t.test(lowGC$meanCAI_a, lowGC$meanCAI_not_a, paired=TRUE)

cat('High GC\n')
t.test(highGC$meanCAI_a, highGC$meanCAI_not_a, paired=TRUE)
cat(sprintf("Mean High GC A: %s +- %s\n", mean(highGC$meanCAI_a), sd(highGC$meanCAI_a)))
cat(sprintf("Mean High GC not A: %s +- %s\n", mean(highGC$meanCAI_not_a), sd(highGC$meanCAI_not_a)))


sink()	


################
# Graphs
################

compress_tiff <- function(filepath){
  library(tiff)
  tiff <- readTIFF(filepath)
  writeTIFF(tiff, filepath, compression="LZW")
}

theme_Publication <- function(base_size=16) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(
      plot.title = element_text(face = "bold" , size = rel(1), hjust =0.5, vjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -3),
      axis.text = element_text(),
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = c(0.9,0.9),
      legend.background = element_rect(fill=NA),
      legend.title=element_blank(),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0")
    ))
}





library(tiff)
library(ggplot2)


# Ecoli codow comparison

# cor <- cor.test(compare_with_codonw$my_cai, compare_with_codonw$codonw_cai, method = "spearman")
# rho1 <- paste(signif(cor$estimate, digits = 6))
# Pval <- paste(signif(cor$p.value, digits =6))
# txt <- bquote(paste(italic(rho), "=", rho1, ", ", italic(P), "=", Pval))
# txt
# txt <- paste("rho == ", "test")

plot <- ggplot(data=compare_with_codonw) +
  geom_point(aes(x=my_cai, y=codonw_cai), col="black") + 
  geom_smooth(aes(x=compare_with_codonw$my_cai, y=compare_with_codonw$codonw_cai),method = "lm", se = FALSE, col="red", lwd=0.5) +
  scale_x_continuous(breaks = seq(0,1,0.1), limits=c(0, 1)) +
  labs(x="CAI values using 20 highly expressed gene indices", y=substitute(paste("CAI values using CodonW ",italic("E. coli"), " indicies"))) +
  theme_Publication()

ggsave('outputs/graphs/8.1_ecoli_codonw_comparison.pdf', plot=plot, dpi=300)
ggsave('outputs/graphs/tiff/8.1_ecoli_codonw_comparison.tiff', plot=plot, dpi=600)
compress_tiff('outputs/graphs/tiff/8.1_ecoli_codonw_comparison.tiff')



data <- data.frame(cai_gc_file$gc3, cai_gc_file$meanCAI)
colnames(data) <- c('gc3', 'cai')
data$group <- ifelse(data$gc3 < 0.2 | data$gc > 0.9, 'extreme', 'regular')
data$group <-as.factor(data$group)

plot <- ggplot(data=data) + 
  geom_point(aes(x=gc3, y=cai, col=group)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits=c(0.0, 1)) +
  scale_colour_manual(name=element_blank(), values=c('2', '1')) +
  labs(x="GC3", y="Genome mean CAI") +
  theme_Publication() + 
  theme(legend.position = "none")

ggsave('outputs/graphs/8.1_GC3_meanCAI.pdf', plot=plot, dpi=300, width=10)
ggsave('outputs/graphs/tiff/8.1_GC3_meanCAI.tiff', plot=plot, dpi=600)
compress_tiff('outputs/graphs/tiff/8.1_GC3_meanCAI.tiff')


print('Outputs in outputs/r_outputs')
print('Graphs in outputs/graphs')



