library(reshape2)
library(PMCMR)




# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))


# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/tiff/"), showWarnings = FALSE)




site_4_ratios <- read.csv('outputs/archaea/ratio_testing/site_4/_site_4_ratios_archaea.csv', head=T)
site_6_ratios <- read.csv('outputs/archaea/ratio_testing/site_6/_site_6_ratios_archaea.csv', head=T)
site_7_ratios <- read.csv('outputs/archaea/ratio_testing/site_7/_site_7_ratios_archaea.csv', head=T)
site_9_ratios <- read.csv('outputs/archaea/ratio_testing/site_9/_site_9_ratios_archaea.csv', head=T)
site_10_ratios <- read.csv('outputs/archaea/ratio_testing/site_10/_site_10_ratios_archaea.csv', head=T)
site_12_ratios <- read.csv('outputs/archaea/ratio_testing/site_12/_site_12_ratios_archaea.csv', head=T)

chitest <- read.csv('outputs/archaea/ratio_testing/site_4/_chisquare_site_4_archaea.csv', head=T)


sink("outputs/r_outputs/14_archaea.txt")

cat("\n===============================================================\n")
cat("Number of genomes with fourth site A > 0.25\n")
cat("===============================================================\n")

genomes_a4 = sum(site_4_ratios$Prop_A4 > 0.25, na.rm=TRUE)
cat(sprintf("%s / %s\n", genomes_a4, length(site_4_ratios$Prop_A4)))
cat(sprintf("%s\n", genomes_a4 / length(site_4_ratios$Prop_A4)))


cat("\n===============================================================\n")
cat("Number of genomes with A4 > 1\n")
cat("===============================================================\n")

genomes_a4_ratio = sum(site_4_ratios$A_Ratio > 1, na.rm=TRUE)
cat(sprintf("%s / %s\n", genomes_a4, length(site_4_ratios$A_Ratio)))
cat(sprintf("%s\n", genomes_a4 / length(site_4_ratios$A_Ratio)))


cat("\n===============================================================\n")
cat("Number of genomes significant enrichment\n")
cat("===============================================================\n")
chitest$padj <- p.adjust(chitest$pval, method="bonferroni")
cat(sprintf('No genomes with significant enrichment (bonferroni correction):  %s/%s\n', sum(ifelse(chitest$padj < 0.05, 1,0)), nrow(site_4_ratios)))




# Are the A ratios at the fourth site correlated with GC3 content?

cat("\n===============================================================\n")
cat("Are the A4 ratios correlated with GC3 content?\n")
cat("===============================================================\n")
shapiro.test(site_4_ratios$GC3)
shapiro.test(site_4_ratios$A_Ratio)

cor.test(site_4_ratios$GC3, site_4_ratios$A_Ratio, method=c("spearman"))


# Are the A4 ratios significantly greater than T4 ratios?

cat("\n=======================================================\n")
cat("Are the A4 ratios significantly greater than T4 ratios?\n")
cat("=======================================================\n")

shapiro.test(site_4_ratios$A_Ratio)
shapiro.test(site_4_ratios$T_Ratio)

wilcox.test(site_4_ratios$A_Ratio, site_4_ratios$T_Ratio, paired=TRUE)


mean_differences <- mean(site_4_ratios$A_Ratio - site_4_ratios$T_Ratio)
cat(sprintf("Mean ratio difference: %s\n\n", mean_differences))

mean_A4_ratios <- mean(site_4_ratios$A_Ratio)
mean_T4_ratios <- mean(site_4_ratios$T_Ratio)

cat(sprintf("Mean A4 ratio: %s\n", mean_A4_ratios))
cat(sprintf("Mean T4 ratio: %s\n", mean_T4_ratios))






# Check for normailty in each of the sites

cat("==============================\n")
cat("Normality checks for each site\n")
cat("==============================\n")

shapiro.test(site_4_ratios$A_Ratio)
shapiro.test(site_7_ratios$A_Ratio)
shapiro.test(site_10_ratios$A_Ratio)


A4ratio <- site_4_ratios$A_Ratio
A6ratio <- site_6_ratios$A_Ratio
A7ratio <- site_7_ratios$A_Ratio
A9ratio <- site_9_ratios$A_Ratio
A10ratio <- site_10_ratios$A_Ratio
A12ratio <- site_12_ratios$A_Ratio


cat("=================================\n")
cat("Variance between synonymous sites\n")
cat("=================================\n")



siteNames <- rep(c("VI", "IX", "XII"), 646)
ratios <- c()

for (i in 1:646) {
  genomeratios <- c(A6ratio[i], A9ratio[i], A12ratio[i])
  ratios <- c(ratios, genomeratios)
}

df <- melt(data.frame(ratios,siteNames), id.vars="siteNames")
kruskal.test(value ~ siteNames, data= df)
posthoc.kruskal.nemenyi.test(x=df$value, g=df$siteNames, method="Chisq")

mean_A6_ratio <- mean(A6ratio)
mean_A9_ratio <- mean(A9ratio)
mean_A12_ratio <- mean(A12ratio)


cat("\nMean Ratios\n")
cat(sprintf("Mean A6 ratio: %s +- %s\n", mean_A6_ratio, sd(A6ratio)))
cat(sprintf("Mean A9 ratio: %s +- %s\n", mean_A9_ratio, sd(A9ratio)))
cat(sprintf("Mean A12 ratio: %s +- %s\n", mean_A12_ratio, sd(A12ratio)))

cat("\n")


cat("=====================================\n")
cat("Variance between non-synonymous sites\n")
cat("=====================================\n")


# Perform kruskall-wallace test
siteNames <- rep(c("IV", "VII", "X"), 646)
ratios <- c()

for (i in 1:646) {
  genomeratios <- c(log(A4ratio[i]), log(A7ratio[i]), log(A10ratio[i]))
  ratios <- c(ratios, genomeratios)
}

df <- melt(data.frame(ratios,siteNames), id.vars="siteNames")
kruskal.test(value ~ siteNames, data= df)
posthoc.kruskal.nemenyi.test(x=df$value, g=df$siteNames, method="Chisq")

mean_A4_ratio <- mean(A4ratio)
mean_A7_ratio <- mean(A7ratio)
mean_A10_ratio <- mean(A10ratio)


cat("\nMean Ratios\n")
cat(sprintf("Mean A4 ratio: %s +- %s\n", mean_A4_ratio, sd(A4ratio)))
cat(sprintf("Mean A7 ratio: %s +- %s\n", mean_A7_ratio, sd(A7ratio)))
cat(sprintf("Mean A10 ratio: %s +- %s\n", mean_A10_ratio, sd(A10ratio)))



sink()


################
# Plots
################
compress_tiff <- function(filepath){
  library(tiff)
  tiff <- readTIFF(filepath)
  writeTIFF(tiff, filepath, compression="LZW")
}

plot <- function(){
  par(mfrow=c(1,1),  xpd=NA)
  boxplot(A4ratio,A7ratio,A10ratio,A6ratio, A9ratio, A12ratio, at=c(1,2,3,5,6,7), boxwex = 0.5, col=c("#a1bfd4","#a1bfd4","#a1bfd4","#cccccc","#cccccc","#cccccc"), cex.axis = 0.9, names =  c("Site 4", "Site 7", "Site 10", "Site 6", "Site 9", "Site 12"), ylab = expression(paste(italic("A"), " Ratio")))
  text( 2, -0.1, "Nonsynonymous sites")
  text( 6, -0.1, "Synonymous sites")
}


pdf("outputs/graphs/14_archaea_site_A_ratios.pdf", width = 30, height=20, pointsize=40)
plot()
dev.off()

tiff("outputs/graphs/tiff/14_archaea_site_A_ratios.tiff", width = 4800, height = 3600, units = "px", res=600)
plot()
dev.off()

compress_tiff('outputs/graphs/tiff/14_archaea_site_A_ratios.tiff')


print('Graphs found in outputs/graphs')
print('Outputs found in outputs/r_outputs/14_archaea.txt')
