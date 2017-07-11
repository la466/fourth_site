
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
site_4_ratios_t4 <- read.csv('outputs/ratio_testing/site_4/_site_4_ratios_t4.csv', head=T)
site_4_ratios_abs <- read.csv('outputs/ratio_testing/site_4/_site_4_abs_A_count.csv', head=T)
site_4_starts <- read.csv('outputs/ratio_testing/site_4/_site_4_start_codon_A_ratios.csv', head=T)
chitest <- read.csv('outputs/ratio_testing/site_4/_chisquare_site_4.csv', head=T)
no_overlaps <- read.csv('outputs/ratio_testing/site_4/_chisquare_site_4_no_overlap.csv', head=T)


################
# Functions
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





################
# Tests
################



sink("outputs/r_outputs/4.1_fourth_site.txt")

cat("\n===============================================================\n")
cat("Number of genomes with fourth site A > 0.25\n")
cat("===============================================================\n")

genomes_a4 <- sum(site_4_ratios$Prop_A4 > 0.25, na.rm=TRUE)
cat(sprintf("%s / %s\n", genomes_a4, length(site_4_ratios$Prop_A4)))

cat("\n===============================================================\n")
cat("Number of genomes with A4 > 1\n")
cat("===============================================================\n")

genomes_a4 <- sum(site_4_ratios$A_Ratio > 1, na.rm=TRUE)
cat(sprintf("%s / %s\n", genomes_a4, length(site_4_ratios$A_Ratio)))
cat(sprintf("%s\n", genomes_a4 / length(site_4_ratios$A_Ratio)))


cat("\n===============================================================\n")
cat("Number of genomes with significant A4 > 1\n")
cat("===============================================================\n")

chitest$padj <- p.adjust(chitest$pval, method="bonferroni")
sum(ifelse(chitest$padj < 0.05, 1, 0))


cat("\n===============================================================\n")
cat("Number of genomes with significant A4 > 1 removing overlaps\n")
cat("===============================================================\n")

no_overlaps$padj <- p.adjust(no_overlaps$pval, method="bonferroni")
sum(ifelse(no_overlaps$padj < 0.05, 1, 0))


cat("\n===============================================================\n")
cat("Max A4 prop\n")
cat("===============================================================\n")

max_a4 <- max(site_4_ratios$Prop_A4, na.rm=TRUE)
cat(sprintf("%s - %s\n", max_a4, site_4_ratios$Genus[site_4_ratios$Prop_A4 == max_a4]))


# Are the A ratios at the fourth site correlated with GC3 content?

cat("\n===============================================================\n")
cat("Are the A4 ratios correlated with GC3 content?\n")
cat("===============================================================\n")
shapiro.test(site_4_ratios$GC3)
shapiro.test(site_4_ratios$A_Ratio)

cor.test(site_4_ratios$GC3, site_4_ratios$A_Ratio, method=c("spearman"))



cat("\n===============================================================\n")
cat("Number of genomes with T4 > 1\n")
cat("===============================================================\n")

genomes_t4 = sum(site_4_ratios$T_Ratio > 1, na.rm=TRUE)
cat(sprintf("%s / %s\n", genomes_t4, length(site_4_ratios$T_Ratio)))
cat(sprintf("%s\n", genomes_t4 / length(site_4_ratios$T_Ratio)))

cat("\n===============================================================\n")
cat("Number of genomes with C4 > 1\n")
cat("===============================================================\n")

genomes_c4 = sum(site_4_ratios$C_Ratio > 1, na.rm=TRUE)
cat(sprintf("%s / %s\n", genomes_c4, length(site_4_ratios$C_Ratio)))
cat(sprintf("%s\n", genomes_c4 / length(site_4_ratios$C_Ratio)))
cat("\n===============================================================\n")
cat("Number of genomes with G4 > 1\n")
cat("===============================================================\n")

genomes_g4 = sum(site_4_ratios$G_Ratio > 1, na.rm=TRUE)
cat(sprintf("%s / %s\n", genomes_g4, length(site_4_ratios$G_Ratio)))
cat(sprintf("%s\n", genomes_g4 / length(site_4_ratios$G_Ratio)))



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





# Start codon use?
cat("\n=======================================================\n")
cat("Start codon A ratios?\n")
cat("=======================================================\n")

cat(sprintf("ATG A4 ratios: %s +- %s\n", mean(site_4_starts$atg_a4_ratio), sd(site_4_starts$atg_a4_ratio)))
cat(sprintf("GTG A4 ratios: %s +- %s\n", mean(site_4_starts$gtg_a4_ratio), sd(site_4_starts$gtg_a4_ratio)))
cat(sprintf("GTG A4 ratios: %s +- %s\n", mean(site_4_starts$ttg_a4_ratio), sd(site_4_starts$ttg_a4_ratio)))


sink()










################
# Plots
################




library(tiff)
library(ggplot2)

# Produce scatter of the A use
lab1 <- expression(paste(italic(A),"-starting second codons", sep=""))
lab2 <- expression(paste(italic(A),"-starting codons", sep=""))

plot <- ggplot(site_4_ratios) +
  geom_point(aes(x=GC3, y=Prop_A4, colour="points1")) +
  geom_smooth(aes(x=GC3, y=Prop_A4,colour="points1"),method = "lm", se = FALSE,  size=0.8) +
  geom_point(aes(x=GC3, y=Prop_A_Codons, colour="points2")) +
  geom_smooth(aes(x=GC3, y=Prop_A_Codons, color="points2"),method = "lm", se = FALSE, size=0.8) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits=c(0.1, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits=c(0.1, 0.7)) +
  labs(x="GC3", y="Genome proportion") +
  scale_colour_manual(name=element_blank(), values=c(points1="blue", points2="black"), labels=c(lab1, lab2)) +
  theme_Publication(base_size=16) +
  theme(legend.position = c(0.25,0.1), legend.text=element_text(size=16), legend.text.align = 0)

ggsave('4.1_A_proportions.pdf', path="outputs/graphs/", plot = plot, dpi=600, width=180, height=180, units="mm")
ggsave('4.1_A_proportions.tiff', path="outputs/graphs/tiff/", plot = plot, dpi=350, width=180, height=180, units="mm")
compress_tiff('outputs/graphs/tiff/4.1_A_proportions.tiff')


# Kernel Density plots
generate_density_plot <- function(site){

  library(tidyr)
  library(dplyr)
  library(ggplot2)
  ratio_file_path <- paste('outputs/ratio_testing/site_', site, '/_site_', site, '_ratios.csv', sep='')

  file <- read.csv(ratio_file_path, head=T)

  aprop <- paste('Prop_A', site, sep='')
  cprop <- paste('Prop_C', site, sep='')
  gprop <- paste('Prop_G', site, sep='')
  tprop <- paste('Prop_T', site, sep='')

  aprops <- file[aprop]
  cprops <- file[cprop]
  gprops <- file[gprop]
  tprops <- file[tprop]

  restrict_columns <- data.frame(aprops, cprops, gprops, tprops)
  colnames(restrict_columns) <- c('A', 'C', 'T', 'G')

  data <- restrict_columns %>% gather()
  colnames(data)<- c("base", 'prop')
  # detach(data)
  # attach(data)
  data.f <- factor(data, levels=c('A', 'C', 'G', 'T'), labels = c("A", "C", "G", "T"))
  cols <- c("#56B4e9", "#e69f00","#c979a7","#009e73")

  data$prop <- as.numeric(data$prop)
  data$base <- as.factor(data$base)

  kplot <- ggplot(data, aes(x=prop)) +
    geom_density(aes(group=base, colour=base), size=1, show.legend=FALSE, lty=1) +
    scale_x_continuous(breaks = seq(0, 1, 0.1), limits=c(0, 1)) +
    labs(x=paste("Proportion of CDSs with nucleotide"), y="Density") +
    stat_density(aes(x=prop, colour=base), geom="line",position="identity") +
    guides(colour = guide_legend(override.aes = list(size=1.2)))+
    scale_colour_manual(values=cols)+
    ggtitle(paste('Site', site)) +
    theme_Publication(base_size=13) +
    theme(legend.position = c(0.85,0.85), legend.title=element_blank())

  return(kplot)
}

plot4 <- generate_density_plot(4)
plot5 <- generate_density_plot(5)
plot6 <- generate_density_plot(6)
plot7 <- generate_density_plot(7)

# save_plot <- function(site, kplot) {
#   ggsave(paste('4.1_site_', site , '_proportion_kernel_density.tiff', sep=""), path="outputs/graphs/tiff/", plot = kplot, dpi=600)
#   tiff <- readTIFF(paste('outputs/graphs/tiff/4.1_site_', site, '_proportion_kernel_density.tiff', sep=""))
#   writeTIFF(tiff, paste('outputs/graphs/tiff/4.1_site_', site, '_proportion_kernel_density.tiff', sep=""), compression="LZW")
# }

# save_plot(4, plot4)
# save_plot(5, plot5)
# save_plot(6, plot6)
# save_plot(7, plot7)


library(gridExtra)
plot <- grid.arrange(plot4, plot5, plot6, plot7, nrow=2)
ggsave('4.1_sites_proportion_kernel_densities.pdf', path="outputs/graphs/", plot = plot, dpi=400)
ggsave('4.1_sites_proportion_kernel_densities.tiff', path="outputs/graphs/tiff/", plot = plot, dpi=600, width=180, height=180, units="mm")
compress_tiff('outputs/graphs/tiff/4.1_sites_proportion_kernel_densities.tiff')


print('Outputs in outputs/r_outputs/4.1_fourth_site.txt')
print('Graphs in outputs/graphs')
