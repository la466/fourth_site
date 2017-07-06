library(ggplot2)
library(reshape2)

# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/tiff/"), showWarnings = FALSE)



magenta <- rep("#b51c81",29)
green <- rep('#26A65B', 29)
blue <- rep('#0072b2', 29)

cbPalette <- c(magenta, green, blue)

################
# File Inputs
################

file <- read.csv("outputs/nucleotide_conservation/_codon_gc_variance.csv", header=TRUE)




colnames(file) <- c("codon", "1", "2","3")


# Melt the data frame to  plot all 3 base positions on the same graph
mdf <- melt(file, id.vars="codon", value.name="value", variable.name="base")

# Rename the columns
colnames(mdf) <- c("codon", "Position", "value")


# Plot the graph of each base position's variance for each codon
plot <- ggplot(data=mdf, aes(x=codon, y=value, group=Position, color=Position)) + geom_point() + geom_line() +
  scale_x_continuous(breaks=c(2:30)) +
  xlab("Codon") + ylab("GC Variance") +
  scale_shape_discrete(name="Experimental\nCondition") 
ggsave("outputs/graphs/5_codon_gc_varaince.pdf", plot = plot, width=10, height=8, dpi=400)
ggsave("outputs/graphs/tiff/5_codon_gc_varaince.tiff", plot = plot, width=10, dpi=400)


compress_tiff <- function(filepath){
  library(tiff)
  tiff <- readTIFF(filepath)
  writeTIFF(tiff, filepath, compression="LZW")
}

compress_tiff('outputs/graphs/tiff/5_codon_gc_varaince.tiff')


print('Graph in outputs/graphs')
