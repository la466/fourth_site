# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))


# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)


library(Cairo)

################
# File Inputs
################

distances <- read.csv("outputs/sd_sequences/_sd_distances.csv", head=T)
sdSequences <- read.csv("outputs/sd_sequences/_sdSequenceAnalysis.csv", head=T)





################
# Tests
################

sink("outputs/r_outputs/13_sd_sequences.txt")


cat("\n===============================================================\n")
cat("Max and min propotion of SD sequences?\n")
cat("===============================================================\n")
max <- max(sdSequences$prop_sd)
min <- min(sdSequences$prop_sd)

cat(sprintf('Max: %s, %s\n', max, sdSequences$genome[sdSequences$prop_sd == max]))
cat(sprintf('Min: %s, %s\n', min, sdSequences$genome[sdSequences$prop_sd == min]))



cat("\n===============================================================\n")
cat("Is the number of SD sequences correlated with GC3?\n")
cat("===============================================================\n")

shapiro.test(sdSequences$gc3)
shapiro.test(sdSequences$prop_sd)

cor.test(sdSequences$gc3, sdSequences$prop_sd, method="spearman")

cat("\n===============================================================\n")
cat("Do CDSs with an SD have more +4A?\n")
cat("===============================================================\n")


shapiro.test(sdSequences$prop_sd_a)
shapiro.test(sdSequences$prop_no_sd_a)
wilcox.test(sdSequences$prop_sd_a, sdSequences$prop_no_sd_a, paired=TRUE)


mean_diff_a_sd_nosd <- mean(sdSequences$prop_sd_a - sdSequences$prop_no_sd_a)
cat(sprintf("Mean A prop SD: %s +- %s\n", mean(sdSequences$prop_sd_a), sd(sdSequences$prop_sd_a)))
cat(sprintf("Mean A prop noSD: %s +- %s\n", mean(sdSequences$prop_no_sd_a), sd(sdSequences$prop_no_sd_a)))
cat(sprintf("Mean difference in A prop SD-noSD: %s\n\n", mean_diff_a_sd_nosd))

length(sdSequences$prop_sd_a)



cat("\n===============================================================\n")
cat("Do CDSs with a weak SD have more +4A?\n")
cat("===============================================================\n")

shapiro.test(sdSequences$prop_weak_sd_a)
shapiro.test(sdSequences$prop_strong_sd_a)
wilcox.test(sdSequences$prop_weak_sd_a, sdSequences$prop_strong_sd_a, paired=TRUE)

mean_diff_a_weak_strong_sd <- mean(sdSequences$prop_weak_sd_a - sdSequences$prop_strong_sd_a)
cat(sprintf("Mean difference in A prop weak-strong: %s\n\n", mean_diff_a_weak_strong_sd))

cat("\n===============================================================\n")
cat("Does the distance of an SD sequence to the initation codon affect +4A?\n")
cat("===============================================================\n")


shapiro.test(sdSequences$near_prop)
shapiro.test(sdSequences$far_prop)
wilcox.test(sdSequences$near_prop, sdSequences$far_prop, paired=TRUE)

mean_diff_a_distance <- mean(sdSequences$near_prop - sdSequences$far_prop)
cat(sprintf("Mean difference in A prop near-far: %s\n\n", mean_diff_a_distance))




sink()

################
# Graphs
################

labels <- c()
for (i in seq(min(distances$pos),max(distances$pos),1)) {
  if (i %% 5 == 0){
    labels <- c(labels,i)
  } else {
    labels <- c(labels,"")
  }
}

label1 <- expression(paste(Delta,"G", degree, " \u2264 ", "-3.4535"))
label2 <- expression(paste("-3.4535 < ", Delta,"G", degree))
label3 <- expression(paste("-8.4 < ",Delta,"G", degree, "  \u2264 -3.4535"))
label4 <- expression(paste(Delta,"G", degree, "  \u2264 -8.4"))
plot_labels <- c(label1, label2, label3, label4)


pdf("outputs/graphs/13_sd_distances.pdf", width=36, height=30, pointsize=40)
par(mfrow=c(3,2))

plot <- barplot(distances$AE000657 / sum(distances$AE000657), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", ylab="Genome proportion", main=substitute(paste(italic('Aquifex'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AE014295 / sum(distances$AE014295), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", main=substitute(paste(italic('Bifidobacterium'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AE005174 / sum(distances$AE005174), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", ylab="Genome proportion", main=substitute(paste(italic('Escherichia'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AL591688 / sum(distances$AL591688), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", main=substitute(paste(italic('Sinorhizobium'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AP012205 / sum(distances$AP012205), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", ylab="Genome proportion", main=substitute(paste(italic('Synechocystis'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AE000512 / sum(distances$AE000512), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", main=substitute(paste(italic('Thermotoga'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

dev.off()

cairo_pdf("outputs/graphs/13_sdSequencesAnalysis.pdf", family="Arial Unicode MS", width = 30,height=20,pointsize = 40)
par(mfrow=c(1,1))
boxplot(sdSequences$prop_sd_a, sdSequences$prop_no_sd_a, sdSequences$prop_weak_sd_a, sdSequences$prop_strong_sd_a,  col = "#a1bfd4", names=plot_labels, ylab=expression(paste("Proportion of genes with +4", italic("A"))), boxwex = 0.5, at=c(1,2,3.5,4.5), cex.axis=0.8, xlab="mRNA - 16S rRNA tail interaction energy (kcal/mol)")
dev.off()





library(tiff)


tiff("outputs/graphs/tiff/13_sd_distances.tiff", width = 8, height=6, units="in", res=400)

par(mfrow=c(3,2))

plot <- barplot(distances$AE000657 / sum(distances$AE000657), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", ylab="Genome proportion", main=substitute(paste(italic('Aquifex'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AE014295 / sum(distances$AE014295), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", main=substitute(paste(italic('Bifidobacterium'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AE005174 / sum(distances$AE005174), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", ylab="Genome proportion", main=substitute(paste(italic('Escherichia'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AL591688 / sum(distances$AL591688), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", main=substitute(paste(italic('Sinorhizobium'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AP012205 / sum(distances$AP012205), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", ylab="Genome proportion", main=substitute(paste(italic('Synechocystis'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

plot <- barplot(distances$AE000512 / sum(distances$AE000512), cex.axis = 1, col = "#a1bfd4", ylim=c(0,0.25), xlab="CDS position", main=substitute(paste(italic('Thermotoga'))))
axis(1,at=plot,label=c(rep("",60)), cex.axis = 0.5)
abline(v=36.7, lty=2, col="#555555")
mtext(side = 1, labels, at=plot, line = 1, cex = 0.5)

dev.off()

## Compress
tiff <- readTIFF("outputs/graphs/tiff/13_sd_distances.tiff")
writeTIFF(tiff, 'outputs/graphs/tiff/13_sd_distances.tiff', compression="LZW")


tiff("outputs/graphs/tiff/13_sdSequencesLowHighSDpropA.tiff", width = 10, height=8, units="in", res=400)
boxplot(sdSequences$prop_sd_a, sdSequences$prop_no_sd_a, sdSequences$prop_weak_sd_a, sdSequences$prop_strong_sd_a,  col = "#a1bfd4", names=plot_labels, ylab=expression(paste("Proportion of genes with +4", italic("A"))), boxwex = 0.5, at=c(1,2,3.5,4.5), cex.axis=0.8, xlab="mRNA - 16S rRNA tail interaction energy (kcal/mol)")
dev.off()

## Compress
tiff <- readTIFF("outputs/graphs/tiff/13_sdSequencesLowHighSDpropA.tiff")
writeTIFF(tiff, 'outputs/graphs/tiff/13_sdSequencesLowHighSDpropA.tiff', compression="LZW")


print('Outputs in outputs/r_outputs')
print('Graphs in outputs/graphs')
