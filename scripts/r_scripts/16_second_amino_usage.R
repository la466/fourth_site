# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))


# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)

################
# File Inputs
################

second_amino_file <- read.csv("outputs/second_amino/_second_amino_ratios.csv", head=TRUE)



################
# Tests
################



sink("outputs/r_outputs/16_second_amino_usage.txt")

mean_serine_a_ratio <- mean(second_amino_file$Sa_ratio)
mean_serine_t_ratio <- mean(second_amino_file$St_ratio)

cat(sprintf("Genomes: %s", length(second_amino_file$Sa_ratio)))

cat(sprintf('Mean serine(A) ratio: %s +- %s\n', mean_serine_a_ratio, sd(second_amino_file$Sa_ratio)))
cat(sprintf('Mean serine(T) ratio: %s +- %s\n', mean_serine_t_ratio, sd(second_amino_file$St_ratio)))

mean_arginine_a_ratio <- mean(second_amino_file$Ra_ratio)
mean_arginine_c_ratio <- mean(second_amino_file$Rc_ratio)


cat(sprintf('Mean arginine(A) ratio: %s +- %s\n', mean_arginine_a_ratio, sd(second_amino_file$Ra_ratio)))
cat(sprintf('Mean arginine(T) ratio: %s +- %s\n', mean_arginine_c_ratio, sd(second_amino_file$Rc_ratio)))


cat("\n=======================================================================================\n")
cat("Is there a significant difference between codon block amino acid usage at the second site?\n")
cat("=======================================================================================\n")

shapiro.test(second_amino_file$Sa_ratio)
shapiro.test(second_amino_file$St_ratio)
wilcox.test(second_amino_file$Sa_ratio, second_amino_file$St_ratio, paired=TRUE)

shapiro.test(second_amino_file$Ra_ratio)
shapiro.test(second_amino_file$Rc_ratio)
wilcox.test(second_amino_file$Ra_ratio, second_amino_file$Rc_ratio, paired=TRUE)



cat("\n=======================================================================================\n")
cat("Compare the usage of amino acids throughout the genome compared with at the second site\n")
cat("=======================================================================================\n")


mean_prop_serine_a_all <- mean(second_amino_file$Sa_prop)
mean_prop_serine_t_all <- mean(second_amino_file$St_prop)


cat(sprintf('Mean prop serine(a) all aminos: %s\n', mean_prop_serine_a_all))
cat(sprintf('Mean prop serine(t) all aminos: %s\n', mean_prop_serine_t_all))
cat(sprintf('serine(A):serine(T) ratio all aminos: %s:%s\n\n', mean_prop_serine_a_all/mean_prop_serine_a_all, mean_prop_serine_t_all/mean_prop_serine_a_all))

mean_prop_serine_a_second <- mean(second_amino_file$Sa_second_prop)
mean_prop_serine_t_second <- mean(second_amino_file$St_second_prop)

cat(sprintf('Mean prop serine(a) second aminos: %s\n', mean_prop_serine_a_second))
cat(sprintf('Mean prop serine(t) second aminos: %s\n', mean_prop_serine_t_second))
cat(sprintf('serine(A):serine(T) ratio all aminos: %s:%s\n\n', mean_prop_serine_a_second/mean_prop_serine_a_second, mean_prop_serine_t_second/mean_prop_serine_a_second))


mean_prop_arginine_a_all <- mean(second_amino_file$Ra_prop)
mean_prop_arginine_t_all <- mean(second_amino_file$Rc_prop)

cat(sprintf('Mean prop arginine(a) all aminos: %s\n', mean_prop_arginine_a_all))
cat(sprintf('Mean prop arginine(t) all aminos: %s\n', mean_prop_arginine_t_all))
cat(sprintf('arginine(A):arginine(T) ratio all aminos: %s:%s\n\n', mean_prop_arginine_a_all/mean_prop_arginine_a_all, mean_prop_arginine_t_all/mean_prop_arginine_a_all))



mean_prop_arginine_a_second <- mean(second_amino_file$Ra_second_prop)
mean_prop_arginine_t_second <- mean(second_amino_file$Rc_second_prop)

cat(sprintf('Mean prop arginine(a) second aminos: %s\n', mean_prop_arginine_a_second))
cat(sprintf('Mean prop arginine(t) second aminos: %s\n', mean_prop_arginine_t_second))
cat(sprintf('arginine(A):arginine(T) ratio all aminos: %s:%s\n\n', mean_prop_arginine_a_second/mean_prop_arginine_a_second, mean_prop_arginine_t_second/mean_prop_arginine_a_second))

sink()

print('Outputs in outputs/r_outputs')
