library(stringi)
library(Cairo)

# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))


# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/tiff/"), showWarnings = FALSE)


cbPalette <- c("#56B4e9", "#e69f00","#CC79A7","#009e73")

################
# File Inputs
################

library(ggplot2)

file <- read.csv("outputs/aod_scores/_aod_second_amino.csv", header=TRUE)

file <- file[file$Base != 'fMa' & file$Base != 'fMc' & file$Base != 'fMg' & file$Base != 'fMt',]
  
as <- c('I', 'M', 'T', 'K', 'Sa', 'Ra')
cs <- c('Rc', 'Lc', 'P', 'H', 'Q')
ts <- c('F', 'Lt', 'St', 'C', 'W', 'Y')
gs <- c('A', 'D', 'E', 'G', 'V')


for (i in seq(1:nrow(file))){
  amino <- file$Base[i]
  if(amino %in% as){
    file$first[i] <- 'A'
  }
  if(amino %in% cs){
    file$first[i] <- 'C'
  }
  if(amino %in% gs){
    file$first[i] <- 'G'
  }
  if(amino %in% ts){
    file$first[i] <- 'T'
  }
}

cols <- c("#56B4e9", "#e69f00", "#c979a7", "#009e73")

file <- file[order(file$first),]
adj <- ifelse(file$low < 0, 1.4,-0.5)

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
      axis.line = element_blank(),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = c(0.9,0.9),
      legend.background = element_rect(fill=NA),
      legend.title=element_text(size = 12),
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


file$Base <- factor(file$Base, levels = file$Base)

adj <- ifelse(file$low < 0, 1.4,-0.4)
plot1 <- ggplot(data=file, aes(x=Base, y=low, fill=first)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=Base), vjust=adj) +
  geom_hline(yintercept = 0, size=0.2) +
  geom_vline(xintercept = 0.4, size=1.3) +
  scale_y_continuous(breaks = seq(-0.03, 0.05, 0.01), limits=c(-0.03, 0.05)) +
  theme_Publication(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.9,0.9),
    plot.title = element_text(face="bold")
  ) +
  labs(y="AOD Score", title=expression(paste("GC", phantom() <= " 0.4419"))) + 
  # scale_colour_manual(name=element_blank(), values=cols) + 
  scale_fill_manual(values=cols, name="First nt")




adj <- ifelse(file$med < 0, 1.4,-0.4)
plot2 <- ggplot(data=file, aes(x=Base, y=med, fill=first)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=Base), vjust=adj) +
  geom_hline(yintercept = 0, size=0.2) +
  geom_vline(xintercept = 0.4, size=1.3) +
  scale_y_continuous(breaks = seq(-0.03, 0.05, 0.01), limits=c(-0.03, 0.05)) +
  theme_Publication(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
  ) +
  labs(y="AOD Score", title=expression(paste("0.4419", phantom() < "GC ", phantom() <= " 0.6091"))) + 
  # scale_colour_manual(name=element_blank(), values=cols) + 
  scale_fill_manual(values=cols, name="First nt")



adj <- ifelse(file$high < 0, 1.4,-0.4)
plot3 <- ggplot(data=file, aes(x=Base, y=high, fill=first)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=Base), vjust=adj) +
  geom_hline(yintercept = 0, size=0.2) +
  geom_vline(xintercept = 0.4, size=1.3) +
  scale_y_continuous(breaks = seq(-0.03, 0.05, 0.01), limits=c(-0.03, 0.05)) +
  theme_Publication(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
  ) +
  labs(y="AOD Score", title=expression(paste("0.6091", phantom() < " GC"))) + 
  # scale_colour_manual(name=element_blank(), values=cols) + 
  scale_fill_manual(values=cols, name="First nt")



library(gridExtra)
plot <- grid.arrange(plot1, plot2, plot3, ncol=1)


compress_tiff <- function(filepath) {
  library(tiff)
  tiff <- readTIFF(filepath)
  writeTIFF(tiff, filepath, compression="LZW")
}

ggsave('outputs/graphs/6_aod.pdf', plot=plot, dpi=600, width=7, height=12)
ggsave('outputs/graphs/tiff/6_aod.tiff', plot = plot, type="cairo", width=180, height=360, dpi=600, units="mm")
compress_tiff('outputs/graphs/tiff/6_aod.tiff')

print('Graph in outputs/graphs')
