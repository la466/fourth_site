# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))


# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)

################
# File Inputs
################

first_stop <- read.csv("outputs/stop_codons/_stopDistanceGenome.csv", header=TRUE)
first_stop_t4 <- read.csv("outputs/stop_codons/_stopDistanceTable4Genome.csv", header=TRUE)

second_stop <- read.csv("outputs/stop_codons/_stopDistanceGenomeSecondStop.csv", header=TRUE)
second_stop_t4 <- read.csv("outputs/stop_codons/_stopDistanceTable4GenomeSecondStop.csv", header=TRUE)

third_stop <- read.csv("outputs/stop_codons/_stopDistanceGenomeThirdStop.csv", header=TRUE)
third_stop_t4 <- read.csv("outputs/stop_codons/_stopDistanceTable4GenomeThirdStop.csv", header=TRUE)


################
# Tests
################


sink("outputs/r_outputs/7_+1_stop_distances.txt")


cat("\n===========================================================\n")
cat("Are the stop distances to the next +1 stop codons different?\n")
cat("===========================================================\n")


cat("\nFirst +1 stop\n")

shapiro.test(first_stop$mean_dist_A)
shapiro.test(first_stop$mean_dist_notA)

wilcox.test(first_stop$mean_dist_A, first_stop$mean_dist_notA, paired = T)

mean_a_first <- mean(first_stop$mean_dist_A)
mean_not_a_first <- mean(first_stop$mean_dist_notA)

cat(sprintf("Mean distance +4A: %s +- %s\n", mean_a_first, sd(first_stop$mean_dist_A)))
cat(sprintf("Mean distance not +4A: %s +- %s\n", mean_not_a_first, sd(first_stop$mean_dist_notA)))
cat(sprintf("Difference: %s\n", mean_a_first-mean_not_a_first))

cat("\nSecond +1 stop\n")

shapiro.test(second_stop$mean_dist_A)
shapiro.test(second_stop$mean_dist_notA)

wilcox.test(second_stop$mean_dist_A, second_stop$mean_dist_notA, paired = T)

mean_a_second <- mean(second_stop$mean_dist_A)
mean_not_a_second <- mean(second_stop$mean_dist_notA)


cat(sprintf("Mean distance +4A: %s +- %s\n", mean_a_second, sd(second_stop$mean_dist_A)))
cat(sprintf("Mean distance not +4A: %s +- %s\n", mean_not_a_second, sd(second_stop$mean_dist_notA)))
cat(sprintf("Difference: %s\n", mean_a_second-mean_not_a_second))


cat("\nThird +1 stop\n")


shapiro.test(third_stop$mean_dist_A)
shapiro.test(third_stop$mean_dist_notA)

wilcox.test(third_stop$mean_dist_A, third_stop$mean_dist_notA, paired = T)

mean_a_third <- mean(third_stop$mean_dist_A)
mean_not_a_third <- mean(third_stop$mean_dist_notA)



cat(sprintf("Mean distance +4A: %s +- %s\n", mean_a_third, sd(third_stop$mean_dist_A)))
cat(sprintf("Mean distance not +4A: %s +- %s\n", mean_not_a_third, sd(third_stop$mean_dist_notA)))
cat(sprintf("Difference: %s\n", mean_a_third-mean_not_a_third))

cat("\n==========\n")
cat("T4 genomes\n")
cat("==========\n")

cat("\nFirst +1 stop_t4\n")

shapiro.test(first_stop_t4$mean_dist_A)
shapiro.test(first_stop_t4$mean_dist_notA)

t.test(first_stop_t4$mean_dist_A, first_stop_t4$mean_dist_notA, paired = T)

mean_a_first <- mean(first_stop_t4$mean_dist_A)
mean_not_a_first <- mean(first_stop_t4$mean_dist_notA)

cat(sprintf("Mean distance +4A: %s\n", mean_a_first))
cat(sprintf("Mean distance not +4A: %s\n", mean_not_a_first))
cat(sprintf("Difference: %s\n", mean_a_first-mean_not_a_first))




sink()


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

library(ggplot2)

plot <- ggplot() +
  geom_point(data = first_stop, aes(x=gc, y=mean_dist)) + 
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.1), limits=c(0.2, 0.8)) +
  labs(x=paste("GC"), y="Mean nucleotide distance to next +1 stop codon") +
  geom_smooth(aes(x=first_stop$gc, y=first_stop$mean_dist),method = "lm", se = FALSE,  size=0.8, col="red") +
  theme_Publication()


ggsave('outputs/graphs/7_mean_distance_+1_stop.pdf', plot=plot, dpi=600)
ggsave('outputs/graphs/tiff/7_mean_distance_+1_stop.tiff', plot = plot, dpi=600)
compress_tiff('outputs/graphs/tiff/7_mean_distance_+1_stop.tiff')





print('Outputs in outputs/r_outputs/7_+1_stop_distances.txt')


