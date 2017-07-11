# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))


# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/tiff/"), showWarnings = FALSE)


file <- read.csv('outputs/ratio_testing/site_4/_site_4_start_codon_A_ratios.csv', head=T)

sink('outputs/r_outputs/4.6_start_codon_a_ratios.txt')

cat('ATG mean A ratio\n')
cat(sprintf('%s +- %s\n\n', mean(file$atg_a4_ratio), sd(file$atg_a4_ratio)))
cat('GTG mean A ratio\n')
cat(sprintf('%s +- %s\n\n', mean(file$gtg_a4_ratio), sd(file$gtg_a4_ratio)))
cat('TTG mean A ratio\n')
cat(sprintf('%s +- %s\n\n', mean(file$ttg_a4_ratio), sd(file$ttg_a4_ratio)))
  
sink()

print('Output in outputs/r_outputs')