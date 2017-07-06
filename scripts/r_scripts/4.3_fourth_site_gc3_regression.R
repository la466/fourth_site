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

################
# Tests
################

sink("outputs/r_outputs/4.3_fourth_site_a_regression_gc3.txt")

cat("\n===============================================================\n")
cat("Regression models of GC3 vs A codon usage\n")
cat("===============================================================\n")


lm_a_all_codons <- lm(site_4_ratios$Prop_A_Codons~site_4_ratios$GC3)
lm_a_second_codon <- lm(site_4_ratios$Prop_A4~site_4_ratios$GC3)

lm_a_all_codons
summary(lm_a_second_codon)

lm_a_second_codon
summary(lm_a_all_codons)

lm_all_coef <- summary(lm_a_all_codons)$coefficients[2,1]
lm_all_se <- summary(lm_a_all_codons)$coefficients[2,2]
lm_second_coef <- summary(lm_a_second_codon)$coefficients[2,1]
lm_second_se <- summary(lm_a_second_codon)$coefficients[2,2]

# Z test: prop A second ≠ prop A codons
z <- (lm_second_coef - lm_all_coef) / sqrt((lm_second_se**2) + (lm_all_se**2))
pvalue2sided <- 2*pnorm(-abs(z))

cat(sprintf('\n\nZ test: is the proportion of A codons vs GC3 = proporiton of A second codons vs GC3?\n\n'))
cat(sprintf("Z score: %s\np-value: %s\n", z, pvalue2sided))

sink()

print('Outputs in outputs/r_outputs/4.3_fourth_site_a_regression_gc3.txt')