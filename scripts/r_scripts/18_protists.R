# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))


# Create directory if not already created
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "graphs/tiff/"), showWarnings = FALSE)


protists <- read.csv('outputs/protists/protists.csv', head=T)
bacteria <- read.csv('outputs/ratio_testing/site_4/_site_4_ratios.csv', head=T)
archaea <- read.csv('outputs/archaea/ratio_testing/site_4/_site_4_ratios_archaea.csv', head=T)


sink('outputs/r_outputs/18_protists.txt')



# Paramecium ratio
cat('Paramecium A4 ratio\n')
para <- protists$A4_ratio[grepl('Paramecium', protists$file)]
cat(sprintf('%s\n\n', para))


# compare with bacteria
cat('\nHow many bacteria with A4 ratio > paramecium\n')
bacteria_greater <- nrow(bacteria[bacteria$A_Ratio > para,])
cat(sprintf('%s/%s, %s%%', bacteria_greater, nrow(bacteria), 100*(bacteria_greater/nrow(bacteria))))
wilcox.test(bacteria$A_Ratio, mu=para)

# compare with archaea
cat('\nHow many archaea with A4 ratio > paramecium\n')
archaea_greater <- nrow(archaea[archaea$A_Ratio > para,])
cat(sprintf('%s/%s, %s%%', archaea_greater, nrow(archaea), 100*(archaea_greater/nrow(archaea))))
wilcox.test(archaea$A_Ratio, mu=para)


# max protists
cat('\nMax protist A4\n')
max_protists <- max(protists$A4_ratio)
cat(sprintf('%s\n\n', max_protists))

# how many bacteria grater than protists
cat('\nHow many bacteria with A4 ratio > max protist\n')
bacteria_greater_max <- nrow(bacteria[bacteria$A_Ratio > max_protists,])
cat(sprintf('%s/%s, %s%%', bacteria_greater_max, nrow(bacteria), 100*(bacteria_greater_max/nrow(bacteria))))


sink()

print('Output in outputs/r_outputs')
