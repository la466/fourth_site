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


sites <- function(organism){

    site_4 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_4/_site_4_ratios.csv', sep=""), head=T)
    site_6 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_6/_site_6_ratios.csv', sep=""), head=T)
    site_7 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_7/_site_7_ratios.csv', sep=""), head=T)
    site_9 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_9/_site_9_ratios.csv', sep=""), head=T)
    site_10 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_10/_site_10_ratios.csv', sep=""), head=T)
    site_12 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_12/_site_12_ratios.csv', sep=""), head=T)
  
    a4 <- site_4$A_Ratio
    a6 <- site_6$A_Ratio
    a7 <- site_7$A_Ratio
    a9 <- site_9$A_Ratio
    a10 <- site_10$A_Ratio
    a12 <- site_12$A_Ratio
  
    t4 <- site_4$T_Ratio
    t6 <- site_6$T_Ratio
    t7 <- site_7$T_Ratio
    t9 <- site_9$T_Ratio
    t10 <- site_10$T_Ratio
    t12 <- site_12$T_Ratio
    
    sink(paste("outputs/r_outputs/15_", organism, "_ratios.txt", sep=""))
    cat('A ratios\n')
    cat(sprintf("4: %s\n", a4))
    cat(sprintf("6: %s\n", a6))
    cat(sprintf("7: %s\n", a7))
    cat(sprintf("9: %s\n", a9))
    cat(sprintf("10: %s\n", a10))
    cat(sprintf("12: %s\n", a12))
    
    cat('\nT ratios\n')
    cat(sprintf("4: %s\n", t4))
    cat(sprintf("6: %s\n", t6))
    cat(sprintf("7: %s\n", t7))
    cat(sprintf("9: %s\n", t9))
    cat(sprintf("10: %s\n", t10))
    cat(sprintf("12: %s\n", t12))
         
    sink()
}

sites('saccharomyces_cerevisiae')
sites('caenorhabditis_elegans')



################
# Plots
################

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

graph <- function(organism){
  print(organism)
  if (organism != 'gene_filtering') {
    

    site_4 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_4/_site_4_ratios.csv', sep=""), head=T)
    site_5 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_5/_site_5_ratios.csv', sep=""), head=T)
    site_6 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_6/_site_6_ratios.csv', sep=""), head=T)
    site_7 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_7/_site_7_ratios.csv', sep=""), head=T)
    site_8 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_8/_site_8_ratios.csv', sep=""), head=T)
    site_9 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_9/_site_9_ratios.csv', sep=""), head=T)
    site_10 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_10/_site_10_ratios.csv', sep=""), head=T)
    site_11 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_11/_site_11_ratios.csv', sep=""), head=T)
    site_12 <- read.csv(paste('outputs/eukaryotes/', organism, '/site_12/_site_12_ratios.csv', sep=""), head=T)
    
    site_4$site <- 4
    site_5$site <- 5
    site_6$site <- 6
    site_7$site <- 7
    site_8$site <- 8
    site_9$site <- 9
    site_10$site <- 10
    site_11$site <- 11
    site_12$site <- 12
  
    data4 <- data.frame(site_4$A_Ratio, site_4$C_Ratio, site_4$G_Ratio, site_4$T_Ratio, site_4$site)
    data5 <- data.frame(site_5$A_Ratio, site_5$C_Ratio, site_5$G_Ratio, site_5$T_Ratio, site_5$site)
    data6 <- data.frame(site_6$A_Ratio, site_6$C_Ratio, site_6$G_Ratio, site_6$T_Ratio, site_6$site)
    data7 <- data.frame(site_7$A_Ratio, site_7$C_Ratio, site_7$G_Ratio, site_7$T_Ratio, site_7$site)
    data8 <- data.frame(site_8$A_Ratio, site_8$C_Ratio, site_8$G_Ratio, site_8$T_Ratio, site_8$site)
    data9 <- data.frame(site_9$A_Ratio, site_9$C_Ratio, site_9$G_Ratio, site_9$T_Ratio, site_9$site)
    data10 <- data.frame(site_10$A_Ratio, site_10$C_Ratio, site_10$G_Ratio, site_10$T_Ratio, site_10$site)
    data11 <- data.frame(site_11$A_Ratio, site_11$C_Ratio, site_11$G_Ratio, site_11$T_Ratio, site_11$site)
    data12 <- data.frame(site_12$A_Ratio, site_12$C_Ratio, site_12$G_Ratio, site_12$T_Ratio, site_12$site)
  
    all <- data.frame(A = numeric(), C = numeric(), G = numeric(), T = numeric(), site=numeric())
    
    
    for (i in 4:12){
      data_name <- paste('data', i, sep="")
      file <- eval(parse(text = paste('site_', i, sep="")))
      data_name <- data.frame(file$A_Ratio, file$C_Ratio, file$G_Ratio, file$T_Ratio, file$site)
      colnames(data_name) <- c('A', 'C', 'G', 'T', 'site')
      all <- rbind(all, data_name)
    }
    
    library(reshape2)
    melt <- melt(all, id.vars = c('site'))
  
    genus <- toupper(substring(organism, 1, 1))
    last <- gsub(".*_","",organism)
    species <- paste(genus, '. ', last, sep="")
    
    
    cols <- c("#56B4e9", "#e69f00","#c979a7","#009e73")
    library(ggplot2)
    ggplot() +
      geom_line(dat=melt, aes(x=site, y=value, col=variable), size=1.2) + 
      theme_Publication() + 
      scale_x_continuous(breaks = seq(4, 12, 1), limits=c(4, 12)) + 
      scale_y_continuous(breaks = seq(0.2, 1.8, 0.1), limits=c(0.2, 1.8)) +
      labs(x="CDS position", y="Enrichment ratio") +
      scale_colour_manual(name=element_blank(), values=cols) +
      theme_Publication(base_size=16) +
      geom_hline(yintercept=1, lty=2) +
      ggtitle(species) + 
      theme(plot.title = element_text(face="italic"))

  }
  

}


genomes <- list.files('outputs/eukaryotes/')

plots <- list()
i <- 0

for (genome in genomes) {
  if (genome != 'gene_filtering'){
    i <- i +1
    plots[[i]] <- graph(genome)
  }
}

library(gridExtra)

all_plots <- do.call(grid.arrange,c(plots, ncol=3))
ggsave('outputs/graphs/15_eukaryote_ratios.pdf', all_plots, height=30, width=20, dpi=300)
ggsave('outputs/graphs/tiff/15_eukaryote_ratios.tiff', all_plots, height=30, width=20, dpi=400)
compress_tiff('outputs/graphs/tiff/15_eukaryote_ratios.tiff')


print('Graphs in outputs/graphs')
print('Outputs in outputs/r_outputs')
