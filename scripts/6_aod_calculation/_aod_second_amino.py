#!/usr/bin/python

# Script number: 			6.1
# File: 					1 of 1
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Calculate average of difference (AOD) score for amino acid in the second position.
#							AOD scores give preference or avoidance in that position.
#							Method from Tang et al. (2010) (http://dx.doi.org/10.1016/j.ygeno.2010.04.001)
# Output file(s):			_aod_second_amino.csv, _aod_all_aminos_all_genomes.csv

import sys, os, csv, re, shutil, numpy
from sys import argv


##########################
# VARIABLES
##########################

# Set the directory containing the sorted and parsed genomes
genomesDir = 'genome_extractions/cds/'
newDir = 'outputs/aod_scores/'

testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE000511"]

# Define the bases
codon_map = {"TTT":"F", "TTC":"F", "TTA":"Lt", "TTG":"Lt",
	   "TCT":"St", "TCC":"St", "TCA":"St", "TCG":"St",
	   "TAT":"Y", "TAC":"Y", "TAA":"taastop", "TAG":"tagstop",
	   "TGT":"C", "TGC":"C", "TGA":"tgastop", "TGG":"W",
	   "CTT":"Lc", "CTC":"Lc", "CTA":"Lc", "CTG":"Lc",
	   "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
	   "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
	   "CGT":"Rc", "CGC":"Rc", "CGA":"Rc", "CGG":"Rc",
	   "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
	   "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
	   "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
	   "AGT":"Sa", "AGC":"Sa", "AGA":"Ra", "AGG":"Ra",
	   "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
	   "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
	   "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
	   "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
	   "ATGstart": "fMa",  "CTGstart": "fMc",  "TTGstart": "fMt",  "GTGstart": "fMg"}


GCrange = ['low', 'med', 'high']
AllAminos = ["fMa", "I", "M", "T", "K", "N", "Sa", "Ra", "fMc", "Rc", "Lc", "P", "H", "Q", "fMt", "F", "Lt", "St", "C", "W", "Y", "fMg", "A", "D", "E", "G", "V"]


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

				# if acc in testFiles:
				genes.append(genus)
				genes.append(row[2:])
				genomes[acc] = genes


		numGenomes = len(genomes)

		return genomes, numGenomes



# Define and open the output files
def newFile():


	aodPath = newDir + "_aod_second_amino.csv"
	aodFile = open(aodPath, "w")


	aodValuesPath = newDir + "_aod_all_aminos_all_genomes.csv"
	aodValuesFile = open(aodValuesPath, "w")

	return aodFile, aodValuesFile



def fileHead(file):

	fileHeadLine = ""
	fileHeadLine += "Base,low,med,high\n"
	file.write(fileHeadLine)



def proportion(number, total):

	if total != 0:
		pc = number / float(total)
	else:
		pc = 0
	return pc





def nucleotideCounts(genomes,acc,numGenomes,fileNumber,aodFile, differences, gc_groupings):

	print ('-' * 10)
	print (acc)
	print ("Genome %d of %d" % (fileNumber, numGenomes))

	genomeDir = genomesDir + acc + ".txt"

	genome = genomes[acc][0]
	loci = genomes[acc][1]

	with open(genomeDir) as myfile:
		line = myfile.read()
		genes = line.split("\n\n")
		genes.pop()

		# Set the gene count for the genome to 0
		geneNumber = 0

		nucleotideCount = 0
		GCcount = 0


		# Set up the dictionaries to hold the nucleotide counts
		# for each position in both the second codon and all codons

		secondCodon = {}
		totalCodons = {}

		for amino in codon_map:
			secondCodon[codon_map[amino]] = 0
			totalCodons[codon_map[amino]] = 0


		# For each gene in the file
		for singleGene in genes:

			# Get the gene locus
			locus = re.findall('locus_tag=(.+?);', singleGene)[0]

			# If the gene locus is in the filtered loci that passed filtering
			if locus in loci:



				# Retrieve the gene sequence from the file
				seq = re.findall('\n(.+)', singleGene)[0]
				seq = seq[0:-3]

				# Add the total number of nuleotides to the total count
				nucleotideCount += len(seq)

				# Add the number of GC nucleotides
				gc = re.subn('[GC]', "0", seq)
				GCcount += gc[1]

				# For each codon in the sequence minus the stop codon
				for i in range(0,len(seq)-3,3):

					# Get the codon sequence
					codon = seq[i:i+3]

					# Get the count for nucleotides in each position of the
					# second codon
					if i == 3:
						secondCodon[codon_map[codon]] += 1

					# Get the count for nucleotides in each position for all
					# codons in the sqeuence
					totalCodons[codon_map[codon]] += 1



		# Calcalate the total number of second codons (number of genes)
		# and the number of codons in these genes
		totalGenes = 0
		totalCodonCount = 0
		for amino in codon_map:
			totalGenes += secondCodon[codon_map[amino]]
			totalCodonCount += totalCodons[codon_map[amino]]


		# Calculate the genome GC content
		GC = proportion(GCcount, nucleotideCount)
		gc_groupings.append(GC)

		# Set up the list to hold the differences between second codon
		# and total codons frequency for the genome
		differences[genome] = {}
		differences[genome][GC] = {}
		for amino in codon_map:
			differences[genome][GC][codon_map[amino]] = {}


		# Calculate difference of nucleotide frequency for all positions
		# between the second codon and all codons
		for amino in codon_map:
			differences[genome][GC][codon_map[amino]] = proportion(secondCodon[codon_map[amino]],totalGenes) - proportion(totalCodons[codon_map[amino]], totalCodonCount)


		return differences


def aodAnalysis(differences, numGenomes, aodFile, aodValuesFile, gc_groupings):

    gc_groupings = numpy.sort(gc_groupings)

    lower_bound = int(len(gc_groupings)/3) -1
    upper_bound = 2*int(len(gc_groupings)/3) -1

    # print(lower_bound, upper_bound)

    lower_gc = gc_groupings[lower_bound]
    upper_gc = gc_groupings[upper_bound]



    # print (lower_gc, upper_gc)

    GClow = 0
    GCmed = 0
    GChigh = 0


    totalDifferences = {}
    for amino in codon_map:
    	totalDifferences[codon_map[amino]] = {}
    	for GCcontent in ['low', 'med', 'high']:
    		totalDifferences[codon_map[amino]][GCcontent] = 0


    headerLine = False

    for genome in differences:
    	for GC in differences[genome]:

    		if headerLine == False:
    			for amino in differences[genome][GC]:
    				fileLine = ",%s" % amino
    				aodValuesFile.write(fileLine)
    			aodValuesFile.write("\n")
    			headerLine = True


    		if GC <= lower_gc:
    			GClow += 1

    			fileLine = "%s,%s" % (genome,GC)
    			aodValuesFile.write(fileLine)


    			for amino in differences[genome][GC]:
    				totalDifferences[amino]['low'] += differences[genome][GC][amino]


    				fileLine = ",%s" % differences[genome][GC][amino]
    				aodValuesFile.write(fileLine)
    				# print genome,GC,amino
    		elif  GC > lower_gc and GC <= upper_gc:
    			GCmed += 1

    			fileLine = "%s,%s" % (genome,GC)
    			aodValuesFile.write(fileLine)


    			for amino in differences[genome][GC]:
    				totalDifferences[amino]['med'] += differences[genome][GC][amino]


    				fileLine = ",%s" % differences[genome][GC][amino]
    				aodValuesFile.write(fileLine)

    				# print genome,GC,amino
    		elif GC > upper_gc:

    			GChigh += 1

    			fileLine = "%s,%s" % (genome,GC)
    			aodValuesFile.write(fileLine)


    			for amino in differences[genome][GC]:
    				totalDifferences[amino]['high'] += differences[genome][GC][amino]


    				fileLine = ",%s" % differences[genome][GC][amino]
    				aodValuesFile.write(fileLine)

    		aodValuesFile.write("\n")



    positionLine = ""



    for amino in AllAminos:

    	 positionLine = "%s" % amino
    	 for GCcontent in GCrange:

			 ## not correct GC boundaries after update, but just refer to low med and high GC genomes

    		 if GCcontent == "low":
    			 positionLine += ",%s" % proportion(totalDifferences[amino][GCcontent],GClow)

    		 elif GCcontent == "med":
    			 positionLine += ",%s" % proportion(totalDifferences[amino][GCcontent],GCmed)

    		 elif GCcontent == "high":
    			 positionLine += ",%s" % proportion(totalDifferences[amino][GCcontent],GChigh)


    	 positionLine += "\n"
    	 # print positionLine

    	 aodFile.write(positionLine)

#################################

def main():

    setupDirectories(newDir)

    # Set the number of files to 0
    fileNumber = 0

    # Return in the information on the genomes
    genomes, numGenomes = getFilteredGenes()


    # Open the new files
    aodFile, aodValuesFile = newFile()

    # Write the heading for the new files
    fileHead(aodFile)

    differences = {}

    gc_groupings = []

    fileNumber = 0
    for acc in genomes:
        # if fileNumber < 6:
		# if acc in testFile:
		fileNumber+=1
		if fileNumber:
			nucleotideCounts(genomes,acc,numGenomes,fileNumber,aodFile, differences, gc_groupings)


    aodValuesFile.write("acc,gc")

    aodAnalysis(differences, numGenomes, aodFile, aodValuesFile, gc_groupings)



    # Close both of the output files
    aodFile.close()
    aodValuesFile.close()

#################################


if __name__ == "__main__":
	main()
