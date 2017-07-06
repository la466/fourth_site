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

file <- read.csv("outputs/ratio_testing/site_comparison/_synonymous_sites.csv", header=TRUE)
# file2 <- read.csv("outputs/ratio_testing/site_comparison/_nonsynonymous_sites.csv", header=TRUE)



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

codons <- c(2:120)


codon_syn = c()
codon_non = c()

for (i in 2:120) {
  codon_syn[i-1]<-cor(file$GC3, file[i+4])
}


library(ggplot2)

plot <- ggplot() +
  geom_point(aes(x=codons, y=codon_syn, colour="points1"), col="black") +
  scale_x_continuous(breaks = seq(0, 120, 10), limits=c(0, 120)) +
  labs(x="Codon", y=expression(paste('GC3 v synonymous site ', italic(A), ' proportion correlation (', rho, ')'))) +
  theme_Publication() + 
  theme(legend.position = "none")



pdf("outputs/graphs/4.4_synonymous_site_gc3_correlation.pdf", width=8)
plot
dev.off()

ggsave('outputs/graphs/tiff/4.4_synonymous_site_gc3_correlation.tiff', plot = plot, dpi=400, width=8, height=6)
compress_tiff('outputs/graphs/tiff/4.4_synonymous_site_gc3_correlation.tiff')


print('Graphs in outputs/graphs')


