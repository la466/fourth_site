# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))


# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)


require(reshape2)
require(PMCMR)

################
# File Inputs
################

leaderDistances <- read.csv("outputs/leader_genes/_leader_distances.csv", head = TRUE)
leaderContent <- read.csv("outputs/leader_genes/_leader_fourth_site.csv", head = TRUE)
leaderDistanceGC <- read.csv("outputs/leader_genes/_leader_distances_GC.csv", header = TRUE)
leaderNoLeader <- read.csv("outputs/leader_genes/_leader_leaderless_a_content.csv", head = T)
genomeAprop <- read.csv("outputs/leader_genes/_genomeAPropLeaderDistances.csv", header = T)
genomeDistance <- read.csv("outputs/leader_genes/_genomeLeaderDistances.csv", header = T)
number_leaders <- read.csv("outputs/leader_genes/_number_leaders.csv", head = TRUE)



################
# Tests
################




sink("outputs/r_outputs/10_leader_genes.txt")

cat("=================================================================\n")
cat("Max / min prop of leaders\n")
cat("================================================================\n")

max <- max(number_leaders$prop)
min <- min(number_leaders$prop)

cat(sprintf('Max prop: %s, %s\n', max, number_leaders$acc[number_leaders$prop == max]))
cat(sprintf('Min prop: %s, %s\n', min, number_leaders$acc[number_leaders$prop == min]))

cat("=================================================================\n")
cat("Does the prop of leader genes corrlate with GC content?\n")
cat("================================================================\n")


# shapiro.test(number_leaders$prop)
# cor.test(number_leaders$gc, number_leaders$prop, method = "spearman")

# shapiro.test(leaderDistanceGC$gc)
# cor.test(leaderDistanceGC$gc, leaderDistanceGC$mean_leader_dist, method = "spearman")


cat("=================================================================\n")
cat("Do leaderless genes have different +4A content?\n")
cat("================================================================\n")


# Determine whether genes with no leader gene have a different 
# +4A content than genes with a potential leader gene

leader_a <- leaderNoLeader$leader_a
no_leader_a <- leaderNoLeader$no_leader_a

shapiro.test(leader_a)
shapiro.test(no_leader_a)

wilcox.test(leader_a, no_leader_a, paired = T)

mean_leader_a <- mean(leader_a)
mean_leaderless_a <- mean(no_leader_a)
mean_diff_a <- mean(leader_a - no_leader_a)

cat(sprintf('Mean +4A prop leader genes: %s +- %s\n', mean_leader_a, sd(leader_a)))
cat(sprintf('Mean +4A prop leaderless genes: %s +- %s\n', mean_leaderless_a, sd(no_leader_a)))
cat(sprintf('Difference: %s\n', mean_diff_a))
cat('\n\n')

cat("=================================================================\n")
cat("Does the distance to the leader gene correlate with +4A content?\n")
cat("================================================================\n")

fourthcontent = leaderContent[1:100,]

# Does leader distance correlate with A content?
cor.test(fourthcontent$distance, fourthcontent$a_content, method="spearman")

cat("=============================================================================\n")
cat("Is +4A content different in genes closer than most common nt distance than further?\n")
cat("============================================================================\n")


# is there a difference in proportion between the fourth site A content 
# in leader genes closer than 13 nucleotides and greater than 13?
fourthcontent = leaderContent[1:100,]

lowerfourthcontent = leaderContent[1:10,]
upperfourthcontent = leaderContent[11:100,]

shapiro.test(lowerfourthcontent$a_content)
shapiro.test(upperfourthcontent$a_content)

t.test(lowerfourthcontent$a_content, upperfourthcontent$a_content)



cat("=============================================================================\n")
cat("Is there a correlation between GC and mean leader distance?\n")
cat("============================================================================\n")


# Is there a correlation between the mean distance between a leader gene
# and the CDS and the GC content?
# shapiro.test(leaderDistanceGC$mean_leader_dist)
# cor.test(leaderDistanceGC$gc3, leaderDistanceGC$mean_leader_dist, method = "spearman")
# 
# shapiro.test(leaderDistanceGC$gc)
# cor.test(leaderDistanceGC$gc, leaderDistanceGC$mean_leader_dist, method = "spearman")
# 
# head(leaderDistanceGC)

sink()





###########
# Graphs
###########


# Leader gene distances
data <- subset(leaderDistances, leaderDistances$acc == "total", select = c(X1:X100))
cols <- c("total")
tdata <- t(data)
colnames(tdata) <- cols
tdata <- data.frame(tdata)

pdf("outputs/graphs/10_leader_gene_distances.pdf", height=24, width=30, pointsize=40)
par(cex.axis=0.8)
par(mfrow=c(1,1))
par(xaxp  = c(0, 200, 2))
barplot(tdata$total,cex = 1, space = 0, xaxt = "n", ylim = c(0,10000), cex.axis = 0.5, col="#a1bfd4", border = "#ffffff", xlab = "Nucleotide distance between the leader sequence stop codon and CDS start codon", ylab="Number of leader genes")
at <- seq(from = 0, to = 100, by = 1)
axis(1,at = c(1:100) - 0.5,labels = c(1:100), cex.axis=0.5)
dev.off()

# Leader gene proportion +4A content with distance
pdf("outputs/graphs/10_leaderGeneDistancesPropA.pdf", height=24, width=30, pointsize=40)
plot(fourthcontent$distance, fourthcontent$a_content, ylim=c(0,0.5), cex=0.6, pch=16, xlab="Nucleotide distance between the leader sequence stop codon and CDS start codon", ylab=expression(paste("Proportion of CDS with +4", italic("A"))), xaxt="na")
axis(1,at = c(1:100) - 0.5,labels = c(1:100), cex.axis=0.5)
dev.off()



library(tiff)
# High res tiffs

tiff("outputs/graphs/tiff/10_leader_gene_distances.tiff", width = 12, height=8, units="in", res=300)
par(cex.axis=0.8)
par(mfrow=c(1,1))
par(xaxp  = c(0, 200, 2))
barplot(tdata$total,cex = 1, space = 0, xaxt = "n", ylim = c(0,10000), cex.axis = 0.5, col="#a1bfd4", border = "#ffffff", xlab = "Nucleotide distance between the leader sequence stop codon and CDS start codon", ylab="Number of leader genes")
at <- seq(from = 0, to = 100, by = 1)
axis(1,at = c(1:100) - 0.5,labels = c(1:100), cex.axis=0.5)
dev.off()

## Compress
tiff <- readTIFF("outputs/graphs/tiff/10_leader_gene_distances.tiff")
writeTIFF(tiff, 'outputs/graphs/tiff/10_leader_gene_distances.tiff', compression="LZW")


tiff("outputs/graphs/tiff/10_leaderGeneDistancesPropA.tiff", width = 12, height=8, units="in", res=300)
plot(fourthcontent$distance, fourthcontent$a_content, ylim=c(0,0.5), cex=0.6, pch=16, xlab="Nucleotide distance between the leader sequence stop codon and CDS start codon", ylab=expression(paste("Proportion of CDS with +4", italic("A"))), xaxt="na")
axis(1,at = c(1:100) - 0.5,labels = c(1:100), cex.axis=0.5)
dev.off()

## Compress
tiff <- readTIFF("outputs/graphs/tiff/10_leaderGeneDistancesPropA.tiff")
writeTIFF(tiff, 'outputs/graphs/tiff/10_leaderGeneDistancesPropA.tiff', compression="LZW")


print('Outputs in outputs/r_outputs')
print('Graphs in outputs/graphs')

