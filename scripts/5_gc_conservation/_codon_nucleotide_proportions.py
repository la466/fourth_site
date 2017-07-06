#!/usr/bin/python

# Script number: 			5.1
# File: 					1 of 3
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv, table4genomes.txt
# Description: 				Caculate the proportions of each nucleotie at each position of codon 2-30.
#							Only use genes that are 30+ codons long
# Output file(s):			_codon_nucleotide_proportions.csv



import sys, os, csv, re, shutil



##########################
# VARIABLES
##########################

# Set the directory containing the sorted and parsed genomes
genomesDir = 'genome_extractions/'

newDir = 'outputs/nucleotide_conservation/'
table4Genomes = 'outputs/gene_filtering/table4genomes.txt'
filteredGenesFile = "outputs/gene_filtering/_filteredGenes.csv"

testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE000511"]

# Define the bases
bases = ["A", "C", "T", "G"]

codonRange = range(2,31)
basePos = range(1,4)


##########################
# FUNCTIONS
##########################

# Set up the new directories for the good files and the bad files
def setupDirectories(direc):

	# If the directory already exists, make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)


def getTable4(file):

	t4genomes = []

	with open(file, "U") as openFile:

		for line in openFile:

			linesplit = line.split(',')
			acc = linesplit[1].rstrip('\n')

			t4genomes.append(acc)

	return t4genomes

# Get all the good genes
def getFilteredGenes():

	goodGenes = {}

	
	with open(filteredGenesFile, 'U') as myfile:
		read = csv.reader(myfile)
		
		genomes = {}

		rowNum = 0
		for row in read:		
			rowNum += 1
			if rowNum > 1:
				
				genus = row[0]
				acc = row[1]

				genes = []

				# if acc in testFile: 
				genes.append(genus)
				genes.append(row[2:])
				genomes[acc] = genes

			
		numGenomes = len(genomes)

		return genomes, numGenomes



# Define and open the output files
def newFile():

	# Output file with genes passing the filter
	nucleotidePath = newDir + "_codon_nucleotide_proportions.csv"
	nucleotideFile = open(nucleotidePath, "w")

	return nucleotideFile



def fileHead(file):

	codonLine = ""
	codonLine += "genus,acc,num_genes,genome_gc"

	for codon in codonRange:
		for pos in basePos:
			for base in bases:
				codonLine += ",codon%s_base%s_%s" % (codon, pos, base)
	
	codonLine += ",codons_nucleotides\n"
	file.write(codonLine)



def proportion(number, total):

	if total != 0:
		pc = number / float(total)
	else:
		pc = 0
	return pc





def secondAminoAnalysis(genomes,acc,numGenomes,fileNumber,nucleotideFile, t4genomes):
	
	print '-' * 10
	print acc
	print "Genome %d of %d" % (fileNumber, numGenomes)

	genomeDir = 'genome_extractions/cds/' + acc + ".txt"

	genome = genomes[acc][0]
	loci = genomes[acc][1]


	if acc not in t4genomes:

		with open(genomeDir) as myfile:
			line = myfile.read()
			genes = line.split("\n\n")
			genes.pop()

			# Set the gene count for the genome to 0
			geneNumber = 0

			nucleotideCount = 0
			GCcount = 0

			codons = {}
			
			for codon in codonRange:
				codons[codon] = {}
				for position in basePos:
					codons[codon][position] = {}
					for base in bases:
						codons[codon][position][base] = 0


			# For each gene in the file
			for singleGene in genes:



				# Get the gene locus
				locus = re.findall('locus_tag=(.+?);', singleGene)[0]

				# If the gene locus is in the filtered loci that passed filtering
				if locus in loci:


					# Retrieve the gene sequence from the file
					seq = re.findall('\n(.+)', singleGene)[0]
					nucleotideCount += len(seq)
					
					gc = re.subn('[GC]', "0", seq[:-3])
					GCcount += gc[1]

					# If the gene contains more than 30 codons
					if len(seq) > max(codonRange)*3:

						# Increase the count of genes
						geneNumber += 1

						# For codons 2-30
						for i in codonRange:

							codonStart = (i*3)-3
							# Get the sequence of the codon
				 			codon = seq[codonStart:codonStart+3]
				 		
				 			for pos in basePos:

				 				codonPos = pos -1
				 				codons[i][pos][codon[codonPos]] += 1


			genomeGC = proportion(GCcount, nucleotideCount)

			genusLine = ""
			genusLine += "%s,%s,%s,%s" % (genome, acc, geneNumber, genomeGC)		 				


			outputLine = ""



			codonsNucleotideCount = 0
			for codon in codons:		
				for position in codons[codon]:
					totalNucleotides = 0
					for base in codons[codon][position]:				
						totalNucleotides += codons[codon][position][base]
						codonsNucleotideCount += codons[codon][position][base]
					
					for base in codons[codon][position]:
						genusLine += ",%s" % (proportion(codons[codon][position][base], totalNucleotides))

			genusLine += ",%s\n" % codonsNucleotideCount
			nucleotideFile.write(genusLine)



		


#################################

def main():

	# Set the number of files to 0
	fileNumber = 0
	
	# Return in the information on the genomes
	genomes, numGenomes = getFilteredGenes()

	setupDirectories(newDir)

	t4genomes = getTable4(table4Genomes)

	# Open the new files
	nucleotideFile = newFile()

	# Write the heading for the new files
	fileHead(nucleotideFile)

	fileNumber = 0
	for acc in genomes:
		fileNumber+=1
		secondAminoAnalysis(genomes,acc,numGenomes,fileNumber,nucleotideFile, t4genomes)


	# Close both of the output files
	nucleotideFile.close()

#################################

# Initiate the filtering
if __name__ == "__main__":
	main()

