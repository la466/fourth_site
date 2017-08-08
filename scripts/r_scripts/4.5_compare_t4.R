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




sink('outputs/r_outputs/4.5_compare_t4.txt')

cat("A4 ratios\n")
cat(sprintf("%s, %f\n", site_4_ratios_t4$Genus, site_4_ratios_t4$A_Ratio))



# Are the ratios significantly greater than those using translation table 4?
cat("\n=======================================================\n")
cat("Are the A4 ratios significantly greater than A4 ratios for table 4 genomes?\n")
cat("=======================================================\n")

shapiro.test(site_4_ratios_t4$A_Ratio)
wilcox.test(site_4_ratios$A_Ratio, site_4_ratios_t4$A_Ratio)


# Are the ratios significantly greater than those using translation table 4?
cat("\n=======================================================\n")
cat("Are A4 ratios for table 4 genomes still lower after removing the synthetic construct?\n")
cat("=======================================================\n")


site_4_ratios_t4_no_synthetic <- site_4_ratios_t4$A_Ratio[site_4_ratios_t4$Genus != 'Synthetic']
wilcox.test(site_4_ratios$A_Ratio, site_4_ratios_t4_no_synthetic)


loess_analysis <- function(site, smoothing_param) {
  
  t11_path <- paste('outputs/ratio_testing/site_', site, '/_site_', site, '_ratios.csv', sep='')
  t4_path <- paste('outputs/ratio_testing/site_', site, '/_site_', site, '_ratios_t4.csv', sep='')
  
  t11 <- read.csv(t11_path, head=T)
  t4 <- read.csv(t4_path, head=T)
  
  reduced_11 <- data.frame(t11$Acc, t11$GC, t11$trans_table, t11$A_Ratio)
  colnames(reduced_11) <- c('acc', 'gc', 'table', 'a_ratio')
  
  # t4$trans_table <- 4
  reduced_4 <- data.frame(t4$Acc, t4$GC, t4$trans_table, t4$A_Ratio)
  colnames(reduced_4) <- c('acc', 'gc', 'table', 'a_ratio')
  
  all <- rbind(reduced_11, reduced_4)
  
  data.lo <- loess(all$a_ratio~all$gc, parametric = F, span=smoothing_param)
  pdf('outputs/graphs/4.5_loess.pdf', height=5)
  plot(all$gc,all$a_ratio, pch=ifelse(all$table==4, 17,16), cex=ifelse(all$table==4, 0.8,0.7), col=ifelse(all$table == 4, 'blue', 'black'), xlab="GC", ylab=bquote(paste(italic('A'), ''['4'], ' enrichment ratio', sep="")))
  j <- order(all$gc)
  lines(all$gc[j], data.lo$fitted[j], col="red", lwd=2)
  legend(0.25,3.2, legend=c(paste("Table 11 (N = ", nrow(t11), ")", sep=""), paste("Table 4 (N = ", nrow(t4), ")", sep="")),col=c('black', 'blue'), pch=c(16,17), cex=0.8,  box.lty=0)
  dev.off()
  
  resids <- data.lo$residuals
  # print(head(data.lo))
  
  cat("\n=======================================================\n")
  cat("Loess 5 t4 genomes\n")
  cat("=======================================================\n")
  krus <- kruskal.test(resids~all$table)
  print(krus)
  
  data <- data.frame(all$gc[j], all$table[j], data.lo$residuals[j])
  colnames(data) <- c('gc', 'table', 'res')
  cat(sprintf('Mean t4 residual: %s\n', mean(data$res[data$table==4])))
  cat(sprintf('Mean t11 residual: %s\n', mean(data$res[data$table==11])))
  
  
  
  t11_path <- paste('outputs/ratio_testing/site_', site, '/_site_', site, '_ratios.csv', sep='')
  t4_path <- paste('outputs/ratio_testing/site_4_ratios_all_t4_genomes.csv', sep='')
  
  t11 <- read.csv(t11_path, head=T)
  t4 <- read.csv(t4_path, head=T)
  
  reduced_11 <- data.frame(t11$Acc, t11$GC, t11$trans_table, t11$A_Ratio)
  colnames(reduced_11) <- c('acc', 'gc', 'table', 'a_ratio')
  
  t4$trans_table <- 4
  reduced_4 <- data.frame(t4$file, t4$gc, t4$trans_table, t4$A4_ratio)
  colnames(reduced_4) <- c('acc', 'gc', 'table', 'a_ratio')
  
  all <- rbind(reduced_11, reduced_4)
  
  data.lo <- loess(all$a_ratio~all$gc, parametric = F, span=smoothing_param)
  resids <- data.lo$residuals
  
  pdf('outputs/graphs/4.5_loess_all_t4_genomes.pdf', height=10)
  par(mfrow=c(2,1))
  plot(all$gc,all$a_ratio, pch=ifelse(all$table==4, 17,16), cex=ifelse(all$table==4, 0.8,0.7), col=ifelse(all$table == 4, 'blue', 'black'), xlab="GC", ylab=bquote(paste(italic('A'), ''['4'], ' enrichment ratio', sep="")))
  j <- order(all$gc)
  lines(all$gc[j], data.lo$fitted[j], col="red", lwd=2)
  legend(0.25,3.2, legend=c(paste("Table 11 (N = ", nrow(t11), ")", sep=""), paste("Table 4 (N = ", nrow(t4), ")", sep="")),col=c('black', 'blue'), pch=c(16,17), cex=0.8,  box.lty=0)
  mtext('A', at =c(0.15), line=2, font=(face=2), cex=1.2)
  
  
  library(vioplot)
  vioplot(resids[all$table==11], resids[all$table==4], col= "gold", names=c(paste("Table 11 (N = ", nrow(t11), ")", sep=""), paste("Table 4 (N = ", nrow(t4), ")", sep="")))
  title(ylab="Loess regression residuals")
  mtext('B', at =c(0.18), line=2, font=(face=2), cex=1.2)
  
  dev.off()
  
  
  cat("\n=======================================================\n")
  cat("Loess all t4 genomes\n")
  cat("=======================================================\n")
  krus <- kruskal.test(resids~all$table)
  print(krus)
  
  data <- data.frame(all$gc[j], all$table[j], data.lo$residuals[j])
  colnames(data) <- c('gc', 'table', 'res')
  cat(sprintf('Mean t4 residual: %s\n', mean(data$res[data$table==4])))
  cat(sprintf('Mean t11 residual: %s\n', mean(data$res[data$table==11])))

  
  
}

loess_analysis(4, 0.75)
sink()


loess_analysis_genera <- function(site, smoothing_param) {
  
  
  
  
  t11_path <- paste('outputs/ratio_testing/site_', site, '/_site_', site, '_ratios.csv', sep='')
  t4_path <- paste('outputs/ratio_testing/site_4_ratios_all_t4_genomes.csv', sep='')
  
  t11 <- read.csv(t11_path, head=T)
  t4 <- read.csv(t4_path, head=T)
  
  t11$genus <- NA
  reduced_11 <- data.frame(t11$Acc, t11$GC, t11$trans_table, t11$A_Ratio, t11$genus)
  colnames(reduced_11) <- c('acc', 'gc', 'table', 'a_ratio', 'genus')
  
  t4$trans_table <- 4
  reduced_4 <- data.frame(t4$acc, t4$gc, t4$trans_table, t4$A4_ratio, t4$genus)
  colnames(reduced_4) <- c('acc', 'gc', 'table', 'a_ratio', 'genus')
  
  all <- rbind(reduced_11, reduced_4)
  
  data.lo <- loess(all$a_ratio~all$gc, parametric = F, span=smoothing_param)
  resids <- data.lo$residuals
  
  
  pdf('outputs/graphs/4.5_loess_all_t4_genomes_genera.pdf', height=5)
  par(mfrow=c(1,1))
  plot(all$gc,all$a_ratio, pch=ifelse(all$table==4, 17, 16), cex=ifelse(all$table==4, 0.8,0.7), col=ifelse(all$table==4, ifelse(all$genus=="Mycoplasma", 'blue', ifelse(all$genus=="Spiroplasma", '#e69f00', ifelse(all$genus=="Mesoplasma","#c979a7", ifelse(all$genus=="Ureaplasma","#009e73",'#56B4e9')))),'black'), xlab="GC", ylab=bquote(paste(italic('A'), ''['4'], ' enrichment ratio', sep="")))
  j <- order(all$gc)
  lines(all$gc[j], data.lo$fitted[j], col="red", lwd=2)
  legend(0.25,3.2, legend=c("Table 11 genomes", "Mycoplasma", "Spiroplasma", 'Mesoplasma', 'Ureaplasma', 'Other'),col=c('black', 'blue', '#e69f00', '#c979a7', '#009e73', '#56B4e9' ), pch=c(16,17,17,17,17,17), cex=0.8,  box.lty=0)
  mtext('C', at =c(0.15), line=2, font=(face=2), cex=1.2)
  dev.off()
  
}

loess_analysis_genera(4, 0.75)

print('Outputs in outputs/r_outputs')
print('Graphs in outputs/graphs')
