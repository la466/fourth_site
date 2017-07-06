library(PMCMR)

closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)

file <- read.csv('outputs/expression/start_codon_CAI_means.csv', head=T)

as <- data.frame(file$atg_cai)
as$nt <- 'A'
colnames(as) <- c('cai', 'nt')

gs <- data.frame(file$gtg_cai)
gs$nt <- 'G'
colnames(gs) <- c('cai', 'nt')

ts <- data.frame(file$ttg_cai)
ts$nt <- 'T'
colnames(ts) <- c('cai', 'nt')

all <- rbind(as, gs, ts)


sink('outputs/r_outputs/8.2_CAI_starts.txt')

cat('Is there a significant difference in mean CAI dependant on start codon?\n\n')
all <- all[complete.cases(all), ]
all$nt <- as.factor(all$nt)
kruskal.test(all$cai~all$nt)

posthoc.kruskal.nemenyi.test(x=all$cai, g=all$nt)

sink()






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

cols <- c("#56B4e9", "#e69f00","#c979a7")
library(ggplot2)
pdf('outputs/graphs/8.2_cai_starts.pdf')
plot <- ggplot(data = all, aes(x=nt, y=cai, fill=nt)) + 
  stat_boxplot(geom='errorbar',coef=2) +
  geom_boxplot(data = all, aes(x=nt, y=cai, fill=nt)) +
  labs(x=paste("Start codon"), y="Genome mean CAI") +
  scale_x_discrete(labels=c('ATG', 'GTG', 'TTG')) +
  scale_fill_manual(values=cols)+
  theme_Publication() + 
  theme(legend.position = "none")
plot
dev.off()

print('Outputs in outputs/r_outputs')
print('Graphs in outputs/graphs')
