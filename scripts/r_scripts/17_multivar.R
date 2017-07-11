# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))



# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)



## Files
multivar_file  <- read.csv('outputs/multivar/_all.csv', head=TRUE)
reduced <- read.csv('outputs/multivar/_reduced.csv', head=TRUE)
genes <- read.csv('outputs/multivar/_gene_level.csv', head=TRUE)

reduced$trans_table <- as.factor(reduced$trans_table)


sink('outputs/r_outputs/17_multivar.txt')

# Linear model just incorporating CAI
cat('Linear model with just CAI\n')
lmcai <- lm(fourth~cai, data=multivar_file)
summary(lmcai)

# Linear model just incorporating SD sequences
cat('Linear model with just SD sequences\n')
lmsd <- lm(fourth~sd, data=multivar_file)
summary(lmsd)

# Linear model using other variables
cat('Linear model with other variables\n')
lm3 <- lm(fourth ~ a_codons + a_content + leader + trans_table, data=reduced)
summary(lm3)



# A model on the gene level with fourth site A a binomial factor
cat('Linear model on gene level\n')
genes$trans_table <- as.factor(genes$trans_table)
lm <- glm(fourth~a_codons + leader*local_a_content*trans_table, data=genes, family=binomial())
summary(lm)

sink()

print('Output in outputs/r_outputs')



