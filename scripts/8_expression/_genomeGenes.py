#!/usr/bin/python

# Script number: 			8.2
# File: 					2 of 5
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 8.1
# Prerequisite file(s):		_filteredGenes.csv, genomesHighExpGenes.txt
# Description: 				Create 2 files for each genome, 1 with all genes, 1 with highly expressed genes
# Output file(s):			accession_no/accession_no.txt, accession_no/accession_no_he.txt


import sys, os, re, shutil, Bio, csv
from sys import argv




##########################
# VARIABLES
##########################


# Set the directory containing the sorted and parsed genomes
genomesDir = 'genome_extractions/cds/'
expression_dir = 'outputs/expression/'

expression_genome_dir = expression_dir + "/genomes/"

# Set the directory with the genome data
goodDir = 'good_genomes/'



testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
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



# Set up a new directory
def setupDirectories(direc):

	# If the directory doesnt exist, make new directory
	if os.path.exists(direc):
		shutil.rmtree(direc)
	os.makedirs(direc)


# Create a list of all the good files
def getGenomes(direc):

	files = []

	# For each file in the downloaded genomes
	for eachfile in os.listdir(direc):

		# if eachfile.startswith(testFile[0]):

			# If the file ends with the known .embl format
			if eachfile.endswith(".embl"):

				# Append the file to the files list
				files.append(eachfile)

	return files, len(files)


# Get a list of the genomes containing the highly expressed genes
def getHighExpGenomes(high_exp_file):

	highExpGenomes = []



	with open(high_exp_file) as myfile:
		for line in myfile:
			reg = re.findall('(.*)(?=\n)', line)
			reg = reg[0] + ".embl"
			highExpGenomes.append(reg)

	return highExpGenomes


# Split the file and return a list of genes
def splitFile(file):
	with open(file) as myfile:
		line = myfile.read()
		genes = line.split("\n\n")
		genes.pop()
		return genes




def parseFile(genome,genes,fileNumber,numFiles):



	# Get the acc number of the genome
	acc = genome

	# Set up a new folder for the genome
	newFolderPath = expression_genome_dir + acc + "/"
	setupDirectories(newFolderPath)

	# Set the path to the new file which is to be written
	cdsPath = newFolderPath + acc + ".txt"

	# Open the new file
	cdsFile = open(cdsPath, "w")



	# Set up the path to the new file containing highly expressed genes
	highExpPath = newFolderPath + acc + "_he.txt"

	# Open the file to write the highly expressed genes to
	highExpFile = open(highExpPath, "w")


	# Create an empty dictionary to contain the locus tags for the genome
	locusTags = []

	# Set the count for possible unknown genes to 0
	unknownCount = 0

	# Set the count for the highly expressed genes to 0
	highCount = 0



	# Print miscellaneous information to the terminal
	print "-" * 10
	print genome
	print "%s of %s" % (fileNumber, numFiles)

	# For each of the genes in the list of genes in the genome
	for singleGene in genes:

		# Extract the gene
		gene = re.findall('gene=(.*?);', singleGene)[0]

		# If the gene isnt blank
		if gene != '':

			# Get the gene
			gene = str(gene).replace("['",'').replace("']",'')

		else:

			# Otherwise increase the unknown count
			unknownCount += 1

			# Get the gene known to unknown
			gene = "unknown%s" % unknownCount

		# Get the sequence for the gene
		seq = re.findall('(?<=\n).*', singleGene)[0]



		# Prepare the lines to write to the file
		# First line containing the gene name and length of gene
		# Second line containing the gene sequence
		fileDesc = ">%s\t\t\t%s residues\n" % (gene,len(seq))
		fileSeq = "%s\n" % seq


		# Write the gene to the file containing all genes
		## if gene not in highExpGenes or gene not in highExpGenesExtra:


		if gene in highExpGenes and highCount < 20:
			highCount += 1
			highExpFile.write(fileDesc+fileSeq)

		cdsFile.write(fileDesc+fileSeq)





	cdsFile.close()
	highExpFile.close()



#################################

def main():

	# Set up the new directory
	setupDirectories(expression_genome_dir)


	# Create the list of good files and number
	goodFiles, numFiles = getGenomes(goodDir)

	fileNumber = 0



	# Get a list of the genomes containing the highly expressed genes
	# that we want to look at
	high_exp_file = "outputs/expression/genomesHighExpGenes.txt"
	highExpGenomes = getHighExpGenomes(high_exp_file)


	# For each genome in the good genomes
	for genome in goodFiles:

		# If the genomes is in the list of genomes where we have highly
		# expressed gene data
		if genome in highExpGenomes:

			fileNumber += 1
			genome = re.findall('(.*?)(?=.embl)', genome)
			genome = genome[0]

			# Get the filepath to the genome text file
			file = genomesDir + genome + ".txt"

			# Split the file into the seperate genes
			genes = splitFile(file)


			parseFile(genome,genes,fileNumber,len(highExpGenomes))


#################################

if __name__ == "__main__":
	main()
