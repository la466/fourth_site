#!/usr/bin/python

# Script number: 			7.2
# File: 					2 of 3
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Calculate the distance to the second out of frame stop codon for genomes comparing A/not A
# Output file(s):			_stopDistanceGenomeSecondStop.csv, _stopDistanceTable4GenomeSecondStop.csv

import sys, os, csv, re, shutil
from sys import argv


filename = argv


# Set the directory containing the sorted and parsed genomes
genomesDir = 'genome_extractions/cds/'

testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE000511"]

# Define the bases
bases = ["A", "C", "T", "G"]

translation_tables = {}
translation_tables["table4"] = ["TAG", "TAA"]
translation_tables["table11"] = ["TAG", "TAA", "TGA"]

stop_directory = "outputs/stop_codons/"





###################################################

def setupDirectories(direc):


	# If the directory already exists, delete and make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)

# Get all the good genes
def getFilteredGenes():

	goodGenes = {}

	filteredGenesFile = "outputs/gene_filtering/_filteredGenes.csv"
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




	stopDistanceDistFile = stop_directory + "_stopDistanceGenomeSecondStop.csv"
	stopDistanceDist = open(stopDistanceDistFile, "w")


	stopDistanceTable4DistributionFile = stop_directory + "_stopDistanceTable4GenomeSecondStop.csv"
	stopDistanceTable4Dist = open(stopDistanceTable4DistributionFile, "w")



	return stopDistanceDist, stopDistanceTable4Dist




def fileHead(file):

	file.write("acc,gc,prop_A,prop_C,prop_T,prop_G,codons_A,codons_C,codons_T,codons_G,A_ratio,C_ratio,T_ratio,G_ratio,mean_dist,mean_dist_A,mean_dist_notA\n")

def fileHeadDist(file):

	file.write("interval,A,notA\n")



def proportion(number, total):

	if total != 0:
		pc = number / float(total)
	else:
		pc = 0
	return pc





def getStopDistance(genomes,acc,numGenomes,fileNumber, stopDistance,stopDistanceTable4):

	print '-' * 10
	print acc
	print "Genome %d of %d" % (fileNumber, numGenomes)

	genomeDir = genomesDir + acc + ".txt"

	genome = genomes[acc][0]
	loci = genomes[acc][1]


	genome_table = ""

	GC = []

	distances = {}
	distances["A"] = []
	distances["notA"] = []

	noStops = {}
	noStops["A"] = 0
	noStops["notA"] = 0

	codons = {}
	secondCodons = {}
	for base in bases:
		codons[base] = 0
		secondCodons[base] = 0






	with open(genomeDir) as myfile:
		line = myfile.read()
		genes = line.split("\n\n")
		genes.pop()

		# Set the gene count for the genome to 0
		geneNumber = 0
		getInfo = {}
		getInfo["trans_table"] = ""



		# For each gene in the file
		for singleGene in genes:

			# Get the gene locus
			locus = re.findall('locus_tag=(.+?);', singleGene)[0]

			# If the gene locus is in the filtered loci that passed filtering
			if locus in loci:


				# Increase the count of genes
				geneNumber += 1

				# Retrieve the gene sequence from the file
				seq = re.findall('\n(.+)', singleGene)[0]

				# Retrieve the gene sequence from the file
				trans_table = "table%s" % re.findall('trans_table=(\d+?);', singleGene)[0]

				# Get the genome translation table to use for
				# printing to the correct file
				if genome_table == "":
					genome_table = trans_table


				# Increase the count of the 4th position
				secondCodons[seq[3]] += 1

				geneLength = len(seq)-3

				# Get the GC content of the gene
				genegc = re.subn('[GC]', "0", seq)[1]
				GC.append(proportion(genegc, geneLength))

				# Get the count of the first base of each codon
				for i in range(0,len(seq)-3,3):
		 			codon = seq[i:i+3]

		 			if codon[0] in bases:
		 				codons[codon[0]] += 1





				# Get the +1 frameshift sequence
				shiftseq = seq[1:-3]


				# Get the length of the frameshift
				framshiftLength = len(shiftseq)


				# Set the distance to the next +1 stop codon to 0
				stopLength = 0

				nextStop = 0

				# For each codon, determine whether the codon is a stop codon
				# If it is, set stopLength to the distance to the codon and
				# update stop to say that a stop codon has been found
				for i in range(3,len(shiftseq)-2,3):
					codon = shiftseq[i:i+3]
					if codon in translation_tables[trans_table]:

						nextStop += 1

						if nextStop == 2:
							stopLength = i

				if stopLength != 0:
					# stopLength = len(shiftseq)

					if seq[3] == "A":
						distances[seq[3]].append(stopLength)
					else:
						distances["notA"].append(stopLength)




		genomeLine = ""
		genomeLine += "%s,%s" % (acc, proportion(sum(GC), len(GC)))



		totalGenes = 0
		for base in secondCodons:
			totalGenes += secondCodons[base]

		totalCodons = 0
		for base in codons:
			totalCodons += codons[base]


		secondProp = {}
		for base in bases:
			secondProp[base] = proportion(secondCodons[base], totalGenes)
			genomeLine += ",%s" % secondProp[base]

		codonProp = {}
		for base in codons:
			codonProp[base] = proportion(codons[base], totalCodons)
			genomeLine += ",%s" % codonProp[base]

		for base in bases:
			genomeLine += ",%s" % proportion(secondProp[base],codonProp[base])

		meanDistance = {}
		meanDistance["all"] = proportion((sum(distances["A"]) + sum(distances["notA"])) , (len(distances["A"]) + len(distances["notA"])))
		meanDistance["A"] = proportion(sum(distances["A"]), len(distances["A"]))
		meanDistance["notA"] = proportion(sum(distances["notA"]), len(distances["notA"]))


		genomeLine += ",%s,%s,%s" % (meanDistance["all"], meanDistance["A"], meanDistance["notA"])
		genomeLine += "\n"

		if genome_table == "table4":
			stopDistanceTable4.write(genomeLine)
		else:
			stopDistance.write(genomeLine)


#################################

def main():

	# Set the number of files to 0
	fileNumber = 0

	setupDirectories(stop_directory)

	# Return in the information on the genomes
	genomes, numGenomes = getFilteredGenes()


	# Open the new files
	stopDistance, stopDistanceTable4 = newFile()

	# Write the heading for the new files
	fileHead(stopDistance)
	fileHead(stopDistanceTable4)


	fileNumber = 0
	for acc in genomes:
		fileNumber+=1
		getStopDistance(genomes,acc,numGenomes,fileNumber,stopDistance,stopDistanceTable4)






	# Close both of the output files
	stopDistance.close()
	stopDistanceTable4.close()

#################################


if __name__ == "__main__":
	main()
