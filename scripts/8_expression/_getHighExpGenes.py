#!/usr/bin/python

# Script number: 			8.1
# File: 					1 of 5
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv, table4genomes.txt
# Description: 				Create a list of the genomes that contain 20 of the list of highly expressed genes
# Output file(s):			_genomesHighExpGenes.txt


import sys, os, csv, re


##########################
# VARIABLES
##########################


# Set the directory containing the sorted and parsed genomes
genomesDir = 'genome_extractions/cds/'
expression_dir = 'outputs/expression/'

# Two lists of test files
testFiles = ["AE000511", "AE006914", "AM711867", "AE005174","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE005174"]


highExpGenes = [
	"rplA", "rplB", "rplC", "rplD", "rplE", "rplF", "rplI", "rplJ", "rplK", "rplL", "rplM", "rplN", "rplO", "rplP", "rplQ", "rplR", "rplS", "rplT", "rplU",
	"rpl1", "rpl2", "rpl3", "rpl4", "rpl5", "rpl6", "rpl9", "rpl10", "rpl11", "rpl12", "rpl13", "rpl14", "rpl15", "rpl16", "rpl17", "rpl18", "rpl19", "rpl20", "rpl21",
	"rpsB", "rpsC", "rpsD", "rpsE", "rpsF", "rpsG", "rpsH", "rpsI", "rpsJ", "rpsK", "rpsL", "rpsM", "rpsN", "rpsO", "rpsP", "rpsQ", "rpsR", "rpsS", "rpsT", "rpsU",
	"rps2", "rps3", "rps4", "rps5", "rps6", "rps7", "rps8", "rps9", "rps10", "rps11", "rps12", "rps13", "rps14", "rps15", "rps16", "rps17", "rps18", "rps19", "rps20", "rps21"
	]




##########################
# FUNCTIONS
##########################


# Get a list of all the good genomes and the number
def getGenomes():

	genomes = []

	for genome in os.listdir(genomesDir):

		genome = os.path.splitext(genome)[0]
		genomes.append(genome)

	return genomes, len(genomes)


# Set up the new directories for the good files and the bad files
def setupDirectories(direc):

	# If the directory already exists, delete and make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)


# Define and open the output files
def newFile():

	filePath = expression_dir + "genomesHighExpGenes.txt"

	file = open(filePath, "w")

	return file


# Create a list of each of the genes
def splitFile(file):
	with open(file) as myfile:
		line = myfile.read()
		genes = line.split("\n\n")
		genes.pop()
		return genes





# Filter the genes using the defined tests
# Pass the genome name, file number, the total number of genomes and new files to be written
def filterGenes(genome,genes,fileNumber,numGenomes, genomesHighExpGenes):

	# Print the genome name and genome number of all genomes to be filtered
	print '-' * 10
	print genome
	print "Genome %d of %d" % (fileNumber, numGenomes)

	transtable = {}

	getInfo = {}
	getInfo["genus"] = ""
	getInfo["acc"] = ""



	highExpGenesCount = 0

	# For each gene in the genome
	for singleGene in genes:



		# Extract the gene from the file
		gene = re.findall('gene=(.*?);', singleGene)[0]

		if gene != '':
			gene = str(gene).replace("['",'').replace("']",'')

		if gene in highExpGenes:
			# print gene
			highExpGenesCount += 1







	highExpCount = False
	if highExpGenesCount >= 20:

		highExpCount = True
		genomesHighExpGenes.write(genome + "\n")


	return highExpCount


def getT4genomes(t4path):

	t4genomes = []

	with open(t4path, "U") as myfile:

		for line in myfile:
			acc = line.strip("\n")
			t4genomes.append(acc)

	return t4genomes

#################################

def main():

	# Set the number of files to 0
	fileNumber = 0

	# Get the list containing the raw files for the good genomes
	genomes, numGenomes = getGenomes()

	# Get a list of table4 genomes
	t4path = "outputs/gene_filtering/table4genomes.txt"
	t4genomes = getT4genomes(t4path)

	# Set up the new directory to contain the information for gene
	# expression
	setupDirectories(expression_dir)


	# Set up a new file to contain a list of genomes which have
	# annotations on the set of highly expressed genes
	genomesHighExpGenes = newFile()


	# Set the counter to 0 for the number of genomes which are to be used
	highExpCount = 0

	for genome in genomes:

		if genome not in t4genomes:

			fileNumber +=1

			# Get the genome text file
			file = genomesDir + genome + ".txt"

			# Split the file into each of the genes
			genes = splitFile(file)

			# Determine whether the genome contains enough of the highly expressed genes
			highCount = filterGenes(genome,genes,fileNumber,numGenomes, genomesHighExpGenes)

			if highCount == True:
				highExpCount +=1

	# Print a final summary of the run
	print "Total genomes: %s" % highExpCount



	genomesHighExpGenes.close()

#################################


if __name__ == "__main__":
	main()
