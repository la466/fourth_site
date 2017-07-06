library(reshape2)
library(PMCMR)



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


site_4_ratios <- read.csv('outputs/ratio_testing/site_4/_site_4_ratios.csv', head=T)
site_6_ratios <- read.csv('outputs/ratio_testing/site_6/_site_6_ratios.csv', head=T)
site_7_ratios <- read.csv('outputs/ratio_testing/site_7/_site_7_ratios.csv', head=T)
site_9_ratios <- read.csv('outputs/ratio_testing/site_9/_site_9_ratios.csv', head=T)
site_10_ratios <- read.csv('outputs/ratio_testing/site_10/_site_10_ratios.csv', head=T)
site_12_ratios <- read.csv('outputs/ratio_testing/site_12/_site_12_ratios.csv', head=T)
site_13_ratios <- read.csv('outputs/ratio_testing/site_13/_site_13_ratios.csv', head=T)
site_15_ratios <- read.csv('outputs/ratio_testing/site_15/_site_15_ratios.csv', head=T)

site_4_ratios_t4 <- read.csv('outputs/ratio_testing/site_4/_site_4_ratios_t4.csv', head=T)
site_7_ratios_t4 <- read.csv('outputs/ratio_testing/site_7/_site_7_ratios_t4.csv', head=T)
site_10_ratios_t4 <- read.csv('outputs/ratio_testing/site_10/_site_10_ratios_t4.csv', head=T)



################
# Test
################


sink("outputs/r_outputs/4.2_site_comparison.txt")


# Check for normailty in each of the sites

cat("==============================\n")
cat("Normality checks for each site\n")
cat("==============================\n")

shapiro.test(site_4_ratios$A_Ratio)
shapiro.test(site_6_ratios$A_Ratio)
shapiro.test(site_7_ratios$A_Ratio)
shapiro.test(site_9_ratios$A_Ratio)
shapiro.test(site_10_ratios$A_Ratio)
shapiro.test(site_12_ratios$A_Ratio)
shapiro.test(site_13_ratios$A_Ratio)
shapiro.test(site_15_ratios$A_Ratio)





A4ratio <- site_4_ratios$A_Ratio
A6ratio <- site_6_ratios$A_Ratio
A7ratio <- site_7_ratios$A_Ratio
A9ratio <- site_9_ratios$A_Ratio
A10ratio <- site_10_ratios$A_Ratio
A12ratio <- site_12_ratios$A_Ratio
A13ratio <- site_13_ratios$A_Ratio
A15ratio <- site_15_ratios$A_Ratio


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


cat("=================================\n")
cat("Extend to 15th site\n")
cat("=================================\n")


# Largest variance / smallest variance < 4, so do Kruskall-Wallace
# test as variables are not normally distributed

siteNames <- rep(c("VI", "IX", "XII", "XV"), 646)
ratios <- c()

for (i in 1:646) {
  genomeratios <- c(A6ratio[i], A9ratio[i], A12ratio[i], A15ratio[i])
  ratios <- c(ratios, genomeratios)
}

df <- melt(data.frame(ratios,siteNames), id.vars="siteNames")
kruskal.test(value ~ siteNames, data= df)
posthoc.kruskal.nemenyi.test(x=df$value, g=df$siteNames, method="Chisq")


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


cat("\n=====================================\n")
cat("Table 4 genomes\n")
cat("=====================================\n")

shapiro.test(site_4_ratios_t4$A_Ratio)
shapiro.test(site_7_ratios_t4$A_Ratio)
shapiro.test(site_10_ratios_t4$A_Ratio)

A4ratio_t4 <- site_4_ratios_t4$A_Ratio
A7ratio_t4 <- site_7_ratios_t4$A_Ratio
A10ratio_t4 <- site_10_ratios_t4$A_Ratio





siteNames <- rep(c("VI", "IX", "XII"), 5)
ratios <- c()

for (i in 1:5) {
  genomeratios <- c(log(A4ratio_t4[i]), log(A7ratio_t4[i]), log(A10ratio_t4[i]))
  ratios <- c(ratios, genomeratios)
}

df <- melt(data.frame(ratios,siteNames), id.vars="siteNames")
kruskal.test(value ~ siteNames, data= df)


cat("\nAre table 4 A4 ratios significantly different from A7 or A10 ratios?\n")
wilcox.test(site_7_ratios$A_Ratio, site_4_ratios_t4$A_Ratio)
wilcox.test(site_10_ratios$A_Ratio, site_4_ratios_t4$A_Ratio)


mean_a4_t4 <- mean(site_4_ratios_t4$A_Ratio)
mean_a7_t4 <- mean(site_7_ratios_t4$A_Ratio)
mean_a10_t4 <- mean(site_10_ratios_t4$A_Ratio)

cat(sprintf('Mean A4 ratio (T4): %s +- %s\n', mean_a4_t4, sd(site_4_ratios_t4$A_Ratio)))
cat(sprintf('Mean A7 ratio (T4): %s +- %s\n', mean_a7_t4, sd(site_7_ratios_t4$A_Ratio)))
cat(sprintf('Mean A10 ratio (T4): %s +- %s\n', mean_a10_t4, sd(site_10_ratios_t4$A_Ratio)))

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
  text( 2, -0.3, "Nonsynonymous sites")
  text( 6, -0.3, "Synonymous sites")
}

pdf("outputs/graphs/4.2_site_A_ratios.pdf", width = 30, height=20, pointsize=40)
plot()
dev.off()

tiff("outputs/graphs/tiff/4.2_site_A_ratios.tiff", width = 4800, height = 3600, units = "px", res=400)
plot()
dev.off()

compress_tiff('outputs/graphs/tiff/4.2_site_A_ratios.tiff')



print('Outputs in outputs/r_outputs')
print('Graphs in outputs/graphs')