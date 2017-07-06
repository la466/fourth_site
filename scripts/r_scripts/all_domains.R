# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))



# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/tiff/"), showWarnings = FALSE)



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

protists <- read.csv('outputs/protists/protists.csv', head=T)
generate_density_plot_protists <- function(site){
  
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  
  
  
  aprops <- protists$A4_ratio
  cprops <- protists$C4_ratio
  gprops <- protists$G4_ratio
  tprops <- protists$T4_ratio
  
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
  
  plot <- ggplot(data, aes(x=prop)) +
    geom_density(aes(group=base, colour=base), size=1, show.legend=FALSE, lty=1) +
    scale_x_continuous(breaks = seq(0, 3.5, 0.5), limits=c(0, 3.5)) +
    scale_y_continuous(breaks = seq(0, 5, 0.5), limits=c(0, 5)) +
    labs(x=paste("Fourth site enrichment ratio"), y="Density") +
    stat_density(aes(x=prop, colour=base), geom="line",position="identity") +
    guides(colour = guide_legend(override.aes = list(size=1.2)))+
    scale_colour_manual(values=cols)+
    ggtitle('Protists') +
    theme_Publication(base_size=14) +
    theme(legend.position = c(0.85,0.85), legend.title=element_blank())
  
  return(plot)
}

genomes <- list.files('outputs/eukaryotes/')
generate_density_plot_eukaryotes <- function(genomes, site){
  
  data <- data.frame(Species = character(), A = double(), C = double(), G = double(), T = double())
  # colnames(data) <- c('species', 'A', 'C', 'G', 'T')
  
  for(genome in genomes) {
    if (genome != 'gene_filtering'){
      file <- read.csv(paste('outputs/eukaryotes/', genome, '/site_4/_site_4_ratios.csv', sep=""), head=T)
      reduced <- data.frame(file$Genus, file$A_Ratio, file$C_Ratio, file$G_Ratio, file$T_Ratio)
      colnames(reduced) <- c('Species', 'A', 'C', 'G', 'T')
      data <- rbind(data, reduced)
    }
  }
  

  
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  
  
  
 
  data <- data %>% gather()
  colnames(data)<- c("base", 'prop')
  # detach(data)
  # attach(data)
  data.f <- factor(data, levels=c('A', 'C', 'G', 'T'), labels = c("A", "C", "G", "T"))
  cols <- c("#56B4e9", "#e69f00","#c979a7","#009e73")

  data$prop <- as.numeric(data$prop)
  data$base <- as.factor(data$base)

  plot <- ggplot(data, aes(x=prop)) +
    geom_density(aes(group=base, colour=base), size=1, show.legend=FALSE, lty=1) +
    scale_x_continuous(breaks = seq(0, 3.5, 0.5), limits=c(0, 3.5)) +
    scale_y_continuous(breaks = seq(0, 5, 0.5), limits=c(0, 5)) +
    labs(x=paste("Fourth site enrichment ratio"), y="Density") +
    stat_density(aes(x=prop, colour=base), geom="line",position="identity") +
    guides(colour = guide_legend(override.aes = list(size=1.2)))+
    scale_colour_manual(values=cols)+
    ggtitle('Selected eukaryotes') +
    theme_Publication(base_size=14) +
    theme(legend.position = c(0.85,0.85), legend.title=element_blank())
  
  return(plot)
}






generate_density_plot <- function(site){
  
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  
  file <- read.csv(paste('outputs/ratio_testing/site_', site, '/_site_', site, '_ratios.csv', sep=""), head=T)
  
  aprops <- file$A_Ratio
  cprops <- file$C_Ratio
  gprops <- file$G_Ratio
  tprops <- file$T_Ratio
  
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
  
  plot <- ggplot(data, aes(x=prop)) +
    geom_density(aes(group=base, colour=base), size=1, show.legend=FALSE, lty=1) +
    scale_x_continuous(breaks = seq(0, 3.5, 0.5), limits=c(0, 3.5)) +
    scale_y_continuous(breaks = seq(0, 5, 0.5), limits=c(0, 5)) +
    labs(x=paste("Fourth site enrichment ratio"), y="Density") +
    stat_density(aes(x=prop, colour=base), geom="line",position="identity") +
    guides(colour = guide_legend(override.aes = list(size=1.2)))+
    scale_colour_manual(values=cols)+
    ggtitle('Bacteria') +
    theme_Publication(base_size=14) +
    theme(legend.position = c(0.85,0.85), legend.title=element_blank())

  return(plot)
  
}
generate_density_plot_archaea <- function(site){
  
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  
  file <- read.csv(paste('outputs/archaea/ratio_testing/site_', site, '/_site_', site, '_ratios_archaea.csv', sep=""), head=T)
  

  aprops <- file$A_Ratio
  cprops <- file$C_Ratio
  gprops <- file$G_Ratio
  tprops <- file$T_Ratio
  
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
  
  plot <- ggplot(data, aes(x=prop)) +
    geom_density(aes(group=base, colour=base), size=1, show.legend=FALSE, lty=1) +
    scale_x_continuous(breaks = seq(0, 3.5, 0.5), limits=c(0, 3.5)) +
    scale_y_continuous(breaks = seq(0, 5, 0.5), limits=c(0, 5)) +
    labs(x=paste("Fourth site enrichment ratio"), y="Density") +
    stat_density(aes(x=prop, colour=base), geom="line",position="identity") +
    guides(colour = guide_legend(override.aes = list(size=1.2)))+
    scale_colour_manual(values=cols)+
    ggtitle('Archaea') +
    theme_Publication(base_size=14) +
    theme(legend.position = c(0.85,0.85), legend.title=element_blank())
  
  return(plot)
  
}


library(gridExtra)
plot <- grid.arrange(generate_density_plot(4), generate_density_plot_archaea(4), generate_density_plot_eukaryotes(genomes,4), generate_density_plot_protists(4))
ggsave('outputs/graphs/ratio_density_all.pdf', plot = plot, dpi=400)
ggsave('outputs/graphs/tiff/ratio_density_all.tiff', plot = plot, dpi=400)
compress_tiff('outputs/graphs/tiff/ratio_density_all.tiff')


print('Graphs in outputs/graphs')
