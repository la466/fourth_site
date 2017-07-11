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

compare_file <- read.csv("outputs/species_compare/_compare_ecoli_shigella.csv", head=TRUE)


################
# Tests
################

# Is the fourth A convervation significnatly different to downstream A content?


sink("outputs/r_outputs/9_ecoli_shigella_compare.txt")

cat("\n=================================================================================\n")
cat("Is fourth site A content conservation significnatly different to downstream site A content?\n")
cat("=================================================================================\n\n")

compareA <- compare_file$A[compare_file$pos > 4]
pos4A <- compare_file$A[compare_file$pos == 4]

shapiro.test(log(compareA))

t.test(log(compareA), mu=log(pos4A) )


sink()



################
# Graphs
################

pos4 <- compare_file[2,]
a4 <- pos4$A

legendtext <- expression(paste("+4", italic("A"), " change"))



compress_tiff <- function(filepath) {
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
      # axis.title = element_text(face = "bold",size = rel(1)),
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
      # legend.position = "right",
      # legend.direction = "vertical",
      # legend.key.size= unit(0.2, "cm"),
      # legend.margin = unit(0.5, "cm"),
      # legend.title = element_text(face="italic"),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0")
      # strip.text = element_text(face="bold")
    ))
}



library(tidyr)
library(dplyr)
library(ggplot2)

gat <- gather(compare_file, key, value, -pos)
cols <- c("#56B4e9", "#e69f00","#c979a7","#009e73")


plot <- ggplot() +
  geom_line(data=gat, aes(x=pos, y=value, colour=key), size=1.5) +
  scale_x_continuous(breaks = seq(1, 31, 3), limits=c(0, 31)) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.05), limits=c(0, 0.4)) +
  labs(x="CDS position", y="Proportion of CDSs with substitution") +
  scale_colour_manual(name=element_blank(), values=cols) +
  geom_hline(aes(yintercept=compare_file$A[compare_file$pos == 4]), lty=2, size=0.5, colour="#555555") +
  annotate("text", label="italic(A)[4]", parse=TRUE, x = 0.01, y = compare_file$A[compare_file$pos == 4]+0.008) +
  theme_Publication(base_size=16)


ggsave('9_ecoli_shigella_comparison.pdf', path="outputs/graphs/", plot = plot, dpi=400)
ggsave('9_ecoli_shigella_comparison.tiff', path="outputs/graphs/tiff/", plot = plot, dpi=600, width=180, height=180, units="mm")
compress_tiff('outputs/graphs/tiff/9_ecoli_shigella_comparison.tiff')




print('Outputs in outputs/r_outputs/9_ecoli_shigella_compare.txt')
print('Graphs in outputs/graphs')
