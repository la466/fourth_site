#!/usr/bin/python

# Script number: 			14.5
# File: 					5 of 7
# Prerequisite script(s):	14.1, 14.2, 14.3, 14.4
# Prerequisite file(s):		
# Description: 				Filter the genes dependant on whether they meet the criteria:
# 								FILTER 1: Gene must be a multiple of 3
# 								FILTER 2: Gene must contain only ACTG
# 								FILTER 3: Gene must used a standard stop codon
# 								FILTER 4: Gene must have no in frame stop codons
#								Filter 5: Gene must use an NTG start codon
# Output file(s):			_filteredGenes.csv, _badGenes.csv, _startCodons.csv, table4genomes.txt


import sys, os, csv, re

print "\nFiltering the genes in each genome\n"

##########################
# VARIABLES
##########################

# Set the directory containing the sorted and parsed genomes
genomesDir = 'archaea_genomes/cds/'

# Two lists of test files
testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AM167904"]


output_directory = "outputs/"
newDir = output_directory + "archaea/gene_filtering/"


# Define the stop codons for the translation tables
translation_tables = {}
translation_tables["table4"] = ["TAG", "TAA"]
translation_tables["table11"] = ["TAG", "TAA", "TGA"]
translation_tables["table25"] = ["TAG", "TAA"]

# Define the bases
bases = ["A", "T", "G", "C"]

# Define the start codons
startCodons = ["ATG", "CTG", "GTG", "TTG"]
startCodonCount = {}

# Define all codons
codons = [
	"AAA", "AAC", "AAG", "AAT",
	"ACA", "ACC", "ACG", "ACT",
	"AGA", "AGC", "AGG", "AGT",
	"ATA", "ATC", "ATG", "ATT",
	"CAA", "CAC", "CAG", "CAT",
	"CCA", "CCC", "CCG", "CCT",
	"CGA", "CGC", "CGG", "CGT",
	"CTA", "CTC", "CTG", "CTT",
	"GAA", "GAC", "GAG", "GAT",
	"GCA", "GCC", "GCG", "GCT",
	"GGA", "GGC", "GGG", "GGT",
	"GTA", "GTC", "GTG", "GTT",
	"TAA", "TAC", "TAG", "TAT",
	"TCA", "TCC", "TCG", "TCT",
	"TGA", "TGC", "TGG", "TGT",
	"TTA", "TTC", "TTG", "TTT"
]

for codon in codons:
	startCodonCount[codon] = 0

##########################
# FUNCTIONS
##########################

# Get a list of all the good genomes and the number
def getGenomes():
	
	genomes = []

	for genome in os.listdir(genomesDir):

		if genome != '.DS_Store':			
			genome = os.path.splitext(genome)[0]
			# if genome in testFiles:
			genomes.append(genome)
					
	return genomes, len(genomes)


# Set up the new directories for the good files and the bad files
def setupDirectories(direc):

	# If the directory doesnt exist, make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)


# Define and open the output files
def newFile():

	# Output file with genes passing the filter
	filteredFile = newDir + "_filteredGenes_archaea.csv"
	filtered = open(filteredFile, "w")

	# Output file containing the genes failing the filter
	badFile = newDir + "_badGenes_archaea.csv"
	bad = open(badFile, "w")

	# Output file containing the genes failing the filter
	startCodonFilePath = newDir + "_startCodons_archaea.csv"
	startCodonFile = open(startCodonFilePath, "w")

	table4genomesPath = newDir + "table4genomes_archaea.txt"
	table4genomesFile = open(table4genomesPath, "w")

	return filtered, bad, startCodonFile, table4genomesFile


# Set up the file heading for the good genes
def fileHeadFiltered(file):	
	filteredHead = "Genus,Acc,Gene\n"
	file.write(filteredHead)


# Set up the file heading for the bad genes
def fileHeadBad(file):	
	badFileHead = "Genus,Genome,Gene,fail: 1=used;2=!x3;3=!actg;4=bas_stop;5=inframe_stop\n"
	file.write(badFileHead)

# Set up the file heading for the bad genes
def fileHeadStarts(file):	
	startHead = "start_codon,number,proportion\n"
	file.write(startHead)


# Create a list of each of the genes
def splitFile(file):
	with open(file) as myfile:
		line = myfile.read()
		genes = line.split("\n\n")
		genes.pop()
		return genes

# Test the length of the gene from the raw sequence
# If the gene length is not a multiple of 3, return True
# Assume the gene will fail the test unless proven otherwise
def geneLen(seq):

	# Set the gene length fail test to True
	lengthFail = True

	# Get the gene length
	geneLength = len(seq)

	# Check for a multiple of 3
	if geneLength % 3 == 0:

		# Set the length fail test to False if the gene passes
		lengthFail = False

	# Return the gene length and outcome of test
	return geneLength,lengthFail


# Test to see whether there are any non ACTG bases in the sequence
# If there is a non ACTG base, return true
# Assume the gene will fail the test unless proven otherwise
def notACTG(seq):

	# Set the gene ACTG fail test to True
	ACTGfail = True

	# Replace all non ACTG characters with 0
	# Function will return the number of replacement
	nonACTG = re.subn('[^actgACTG]', '0',seq)

	# If the number of replacements is 0
	if nonACTG[1] == 0:

		# Set the ACTG fail test to False if the gene passes
		ACTGfail = False

	# Return the result of the test and number of ACTG replacements	
	return ACTGfail, nonACTG[1]


# Test to see whether the gene uses a non standard stop codon
# Check the transaltion table to determine the stop codons
# If the gene uses a non standard codon, return True
# Assume the gene will fail the test unless proven otherwise
def nonStop(trans_table,seq):

	# Set the non standard codon test to True
	notStop = True

	# Set the translation table of the genome
	table = "table"+trans_table

	# Retrieve the stop codon
	stop = seq[-3:]

	# If the stop os a regular stop codon for its translation table
	if stop in translation_tables[table]:

		# Set the stop codon fail test to False
		notStop = False

	# Return the result of the test and the stop codon
	return notStop,stop


# Test to see whether the gene has an in frame stop codon
# Check the transaltion table to determine the stop codons
# If the gene has an in frame stop, return True
def inFrameStop(trans_table,seq):

	# Set the in frame stops to False
	inStop = False

	# Set the counter to 0
	inframestop = 0

	# Set the translation table for the genome
	table = "table"+trans_table
	


	#For each codon in the gene sequence, minus the last codon
	for i in range(0,len(seq)-3,3):

		# Retrieve the codon
		codon = seq[i:i+3]

		# If the codon is in the translation tables stop codons
		if codon in translation_tables[table]:

			# Increase the in frame stop counter
			inframestop += 1

	# If the number of in frame stops is greater than 0, set the in frame stop test to true		
	if inframestop > 0:
		inStop = True

	# Return the result of the test
	return inStop


# Check to see whether the gene uses a standard start codon
# If the gene has aa standard start codon, return True
def nonStart(seq):

	# Assume the start codon is non standard
	nonStartCodon = True

	startCodon = seq[0:3]

	# If it is a standard start, return true
	if startCodon in startCodons:
		nonStartCodon = False

	return nonStartCodon


# Determine the proportion given a number and the total number
def proportion(number, total):
	prop = number / float(total)
	return prop



# Filter the genes using the defined tests
# Pass the genome name, file number, the total number of genomes and new files to be written
def filterGenes(genome,genes,fileNumber,numGenomes,filtered,bad, table4genomes):
	
	# Print the genome name and genome number of all genomes to be filtered
	print '-' * 10
	print genome
	print "Genome %d of %d" % (fileNumber, numGenomes)

	transtable = {}

	getInfo = {}
	getInfo["genus"] = ""
	getInfo["acc"] = ""

	filteredGenes = ''

	goodGenes = []
	badGenes = []

	trans_table4 = False

	# For each gene in the genome
	for singleGene in genes:

		

		# Retrieve the genus from the file
		genus = re.findall('^>(.+?);',singleGene)[0]
		getInfo["genus"] = genus


		# Extract the gene name from the file
		acc = re.findall('acc=(.+?);', singleGene)[0]
		getInfo["acc"] = acc
		
		# Extract the gene name from the file
		gene = re.findall('locus_tag=(.+?);', singleGene)[0]


		# # Retrieve the gene sequence from the file
		seq = re.findall('\n(.+)', singleGene)[0]

		
		# If the gene has not already been filtered
		if gene in goodGenes or gene in badGenes:
			print "gene already used: %s" % gene
			geneUsed = "%s,%s,%s,1\n" % (genus, genome, gene)
			bad.write(geneUsed)
			badGenes.append(gene)
		else:
			# FILTER 1
			# Check to see whether gene is multiple of 3
			geneLength, lengthFail = geneLen(seq)

			# If the gene fails the test, write the genus, gene and error code (2) to bad file
			# If the gene passes, continue
			if lengthFail:
				failMutlipleThree = "%s,%s,%s,2\n" % (genus, genome, gene)
				bad.write(failMutlipleThree)
				badGenes.append(gene)
			else:

				# FILTER 2
				# Check for non ACTG bases in sequence
				ACTGfail, nonACTG = notACTG(seq)

				# If the gene fails the test, write the genus, gene and error code (3) to bad file
				# If the gene passes, continue
				if ACTGfail:
					failACTG = "%s,%s,%s,3\n" % (genus, genome, gene)
					bad.write(failACTG)
					badGenes.append(gene)
				else:

					# FILTER 3
					# Check for non standard stop codon

					# Retrive the translation table for the gene
					trans_table = re.findall('trans_table=\d+',singleGene)[0]
					trans_table = re.findall('\d+', trans_table)[0]
					transtable["table"] = trans_table

					if trans_table == '4' and trans_table4 == False:
						
						tableLine = "%s,%s\n" % (getInfo["genus"], getInfo["acc"]) 
						table4genomes.write(tableLine)
						trans_table4 = True
					
					notStop,stop = nonStop(trans_table,seq)

					# If the gene fails the test, write the genus, gene and error code (4) to bad file
					# If the gene passes, continue
					if notStop:
						failStopCodon = "%s,%s,%s,4\n" % (genus, genome, gene)
						bad.write(failStopCodon)
						badGenes.append(gene)
					else:

						# FILTER 4
						# Check for an inframe stop
						inStop = inFrameStop(trans_table,seq)

						# If the gene fails the test, write the genus, gene and error code (5) to bad file
						# If the gene passes, continue
						if inStop:
							failInStop = "%s,%s,%s,5\n" % (genus, genome, gene)
							bad.write(failInStop)
							badGenes.append(gene)
						else:

							# Get the start codon
							startCodon = seq[0:3]

							# Increase the count for the start codons
							if startCodon in startCodonCount:
								startCodonCount[startCodon] += 1
							else:
								startCodonCount[startCodon] = 1

							# FILTER 5
							# Check for a standard start codon
							nonRegularStart = nonStart(seq)

							# If the gene fails the test, write the genus, gene and error code (6) to bad file
							# If the gene passes, continue
							if nonRegularStart:
								badStart = "%s,%s,%s,6\n" % (genus, genome, gene)
								bad.write(badStart)
								badGenes.append(gene)
							else:

								# The gene has passed the filtering
								passGene = ",%s" % gene
								filteredGenes += passGene
								goodGenes.append(gene)


	filteredLine = getInfo["genus"] + "," + getInfo["acc"] + filteredGenes + "\n"
	filtered.write(filteredLine)		


#################################

def main():

	setupDirectories(output_directory)
	setupDirectories(newDir)

	# Set the number of files to 0
	fileNumber = 0
	
	# Return in the information on the genomes
	genomes, numGenomes = getGenomes()
	
	# Open the new files
	filtered,bad,startCodonFile, table4genomes = newFile()

	# Write the heading for the new files
	fileHeadFiltered(filtered)
	fileHeadBad(bad)
	fileHeadStarts(startCodonFile)

	# Filter the genes
	for genome in genomes:
		fileNumber +=1
		file = genomesDir + genome + ".txt"
		genes = splitFile(file)
		filterGenes(genome,genes,fileNumber,numGenomes,filtered,bad, table4genomes)

	# Calculate the number of genes analysed
	startGenes = 0
	for start in startCodonCount:
		startGenes += startCodonCount[start]
	
	# Write the usage of each start codon to the file
	for start in startCodonCount:
		startLine = "%s,%s,%s\n" % (start, startCodonCount[start], proportion(startCodonCount[start], startGenes))
		startCodonFile.write(startLine)


	# Calculate the sum of the proportions
	proportionSum = 0
	for start in startCodonCount:
		proportionSum += proportion(startCodonCount[start], startGenes)

	startsTotal = "Total,%s,%s" % (startGenes, proportionSum)
	startCodonFile.write(startsTotal)

	# Close the output files
	filtered.close()
	bad.close()
	startCodonFile.close()
	table4genomes.close()

#################################


if __name__ == "__main__":
	main()
