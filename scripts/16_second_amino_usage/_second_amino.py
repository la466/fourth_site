#!/usr/bin/python

# Script number: 			16.2
# File: 					2 of 2
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Get the ratios of amino acid use
# Output file(s):			_second_amino_ratios.csv, second_amino_ratios_t4.csv

import sys, os, csv, re
from sys import argv
import rpy2.robjects as robjects




############
# Variables
############


# Set the directory containing the sorted and parsed genomes
genomesDir = 'bacterial_genomes/codons/genomes/'
output_dir = "outputs/second_amino/"
filteredGenesFile = "outputs/gene_filtering/_filteredGenes.csv"
cds_dir = "genome_extractions/cds/"


testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["CP000925"]

# Define the bases
bases = ["A", "C", "T", "G"]

translation_tables = {}
translation_tables["table4"] = ["TAG", "TAA"]
translation_tables["table11"] = ["TAG", "TAA", "TGA"]
translation_tables["table25"] = ["TAG", "TAA"]

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

codon_map_table4 = {"TTT":"F", "TTC":"F", "TTA":"Lt", "TTG":"Lt",
       "TCT":"St", "TCC":"St", "TCA":"St", "TCG":"St",
       "TAT":"Y", "TAC":"Y", "TAA":"taastop", "TAG":"tagstop",
       "TGT":"C", "TGC":"C", "TGA":"W", "TGG":"W",
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



AstartAminos = ["fMa", "I", "M", "T", "K", "Sa", "Ra", "N"]
CstartAminos = ["fMc", "Lc", "Rc", "P", "H", "Q"]
TstartAminos = ["fMt", "F", "Lt", "St", "C", "W", "Y"]
GstartAminos = ["fMg", "A", "D", "E", "G", "V"]
AllAminos = ["fMa", "I", "M", "T", "K", "N", "Sa", "Ra", "fMc", "Rc", "Lc", "P", "H", "Q", "fMt", "F", "Lt", "St", "C", "W", "Y", "fMg", "A", "D", "E", "G", "V"]




###########
# Functions
###########


# Setup new directories
def setupDirectory(directory):

	if not os.path.exists(directory):
		os.makedirs(directory)
		print "Created the new directory: %s" % directory

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






def fileHead(file):

	file.write("genus,acc,trans_table,num_genes,gc3")

	# write the header for the base information
	for base in bases:
		basesummary = ",total_second_%s,total_%s,%s_ratio" % (base,base,base)
		file.write(basesummary)

	# write the header for the amino acids:
	for amino in AllAminos:
		aminoHead = ",%s_second_obs,%s_second_prop,%s_total,%s_prop,%s_ratio" %(amino,amino,amino,amino,amino)
		file.write(aminoHead)

	file.write(",total_aminos\n")


def proportion(number, total):

	if total != 0:
		pc = number / float(total)
	else:
		pc = 0
	#pc = float("%.8f" % pc)
	return pc





def secondAminoAnalysis(genomes,acc,numGenomes,fileNumber,secondAmino,secondAminoTable4):

	print '-' * 10
	print acc
	print "Genome %d of %d" % (fileNumber, numGenomes)

	genomeDir = cds_dir + acc + ".txt"

	genome = genomes[acc][0]
	loci = genomes[acc][1]

	with open(genomeDir) as myfile:
		line = myfile.read()
		genes = line.split("\n\n")
		genes.pop()

		# Set the gene count for the genome to 0
		geneNumber = 0
		getInfo = {}
		getInfo["trans_table"] = ""


		# Dictionary containing the counts for second amino acid counts
		aminoPosTwoCount = {}
		aminoCount = {}

		# For each of the codons in the codon_map
		for amino in AllAminos:

			# Set the number of amino acids they code for to 0
			aminoPosTwoCount[amino] = 0
			aminoCount[amino] = 0


		GC3 = 0

		totalAStartingPosTwo = 0
		totalPosTwo = 0

		totalStarting = {}
		totalStartingPosTwo = {}

		for base in bases:
			totalStarting[base] = 0
			totalStartingPosTwo[base] = 0



		# For each gene in the file
		for singleGene in genes:

			# Get the gene locus
			locus = re.findall('locus_tag=(.+?);', singleGene)[0]

			# print locus

			# If the gene locus is in the filtere loci that passed filtering
			if locus in loci:

				# Increase the count of genes
				geneNumber += 1

				# Retrieve the gene sequence from the file
				seq = re.findall('\n(.+)', singleGene)[0]

				# Retrieve the gene sequence from the file
				getInfo["trans_table"] = re.findall('trans_table=(\d+?);', singleGene)[0]

				# Set the amino acid string to blank
				aminos = ""

				# Retrieve the second codon
		 		codon2 = seq[3:6]

		 		# Increase the count for the codon starts


	 			# Information for the second codon first base usage
	 			if getInfo["trans_table"] == "4":

	 				# Increase the position two count for the amino
		 			aminoPosTwoCount[codon_map_table4[codon2]] += 1

		 			# Increase the count for second codons
	 				# totalAStartingPosTwo += 1
	 				totalStartingPosTwo[codon2[0]] += 1
		 		else:

		 			# Increase the position two count for the amino
					aminoPosTwoCount[codon_map[codon2]] += 1

					# Increase the count for second codons
	 				# totalAStartingPosTwo += 1
		 			totalStartingPosTwo[codon2[0]] += 1



				# For each codon in the gene
				for i in range(0,len(seq)-3,3):

					# Get the sequence of the codon
		 			codon = seq[i:i+3]

		 			if codon[0] in bases:
 						totalStarting[codon[0]] += 1


		 			# Add to GC3
		 			if codon[2] == "G" or codon[2] == "C":
		 				GC3 += 1

		 			# For the first codon
		 			if i == 0:

		 				startcodon = codon + "start"
		 				aminoCount[codon_map[startcodon]] += 1

		 			# For the other codons
		 			else:

		 				# If the translation table is 4
		 				if getInfo["trans_table"] == "4":

		 					# Increase the overall count of the amino
		 					aminoCount[codon_map_table4[codon]] += 1

	 					# For the other translation tables
		 				else:

		 					# Increase the overall count of the amino
							aminoCount[codon_map[codon]] += 1







		# Calculate the number of second codons
		# which should be equal to the number of genes
		# Calculate the total number of amino acids
	 	secondAminoCount = 0
	 	totalCount = 0
	 	for amino in AllAminos:
	 		secondAminoCount += aminoPosTwoCount[amino]
	 		totalCount += aminoCount[amino]



		# Calculate the ratio for the each base in position 4
		# Calculate the ratio of codons starting with each base
		baseFourRatio = {}
		allBaseRatio = {}
	 	for base in bases:
	 		baseFourRatio[base] = proportion(totalStartingPosTwo[base], secondAminoCount)
	 		allBaseRatio[base] = proportion(totalStarting[base], totalCount)

	 	# Calculate the ratio of bases in position 4 to their usage in all codons
	 	baseRatio = {}
	 	for base in bases:
	 		baseRatio[base] = proportion(baseFourRatio[base], allBaseRatio[base])



	 	# Calculate the ratio of second amino acid to total second amino acids
	 	# Calculate the ratio of amino acid to total amino acids
	 	count = {}
	 	for base in bases:
	 		count[base] = 0

	 	secondRatio = {}
	 	totalRatio = {}
	 	for amino in AllAminos:
	 		secondRatio[amino] = proportion(aminoPosTwoCount[amino], secondAminoCount)
	 		totalRatio[amino] = proportion(aminoCount[amino], totalCount)



		# Calculate the second amino acid ratio in relation to whole amino acids

	 	aminoRatio = {}
	 	for amino in AllAminos:
	 		aminoRatio[amino] = proportion(secondRatio[amino], totalRatio[amino])

		# Miscellaneous  information on the genes
		summaryMisc = "%s,%s,%s,%s,%s" % (genome,acc,getInfo["trans_table"], geneNumber, proportion(GC3, totalCount))

		# Write the total second codons starting with each base, total codons starting with each base
		# and the ratio for ach base
		summaryMiscBase = ''
		for base in bases:
			summaryMiscBase += ",%s,%s,%s" % (totalStartingPosTwo[base], totalStarting[base], baseRatio[base])

		# For each amino acid, write the number as codon 2, proportion of codon 2 usage, total useage of the
		# codon, proportion of that amino usage, ratio between codon 2 usage and overall usage
		summaryAmino = ""
		for amino in AllAminos:
			summaryAmino += ",%s,%s,%s,%s,%s" % (aminoPosTwoCount[amino],secondRatio[amino],aminoCount[amino],totalRatio[amino],aminoRatio[amino])

		summaryEnd = ",%s\n" % totalCount

		if getInfo["trans_table"] == "11":
			secondAmino.write(summaryMisc + summaryMiscBase + summaryAmino + summaryEnd)
		else:
			secondAminoTable4.write(summaryMisc + summaryMiscBase + summaryAmino + summaryEnd)








#################################

def main():

	# Set the number of files to 0
	fileNumber = 0

	# Return in the information on the genomes
	genomes, numGenomes = getFilteredGenes()


	setupDirectory(output_dir)

	# Open the new files
	secondAminoFile = output_dir + "_second_amino_ratios.csv"
	secondAmino = open(secondAminoFile, "w")

	# Output file containing the genes failing the filter
	secondAminoTable4File = output_dir + "_second_amino_ratios_t4.csv"
	secondAminoTable4 = open(secondAminoTable4File, "w")

	# Write the heading for the new files
	fileHead(secondAmino)
	fileHead(secondAminoTable4)

	fileNumber = 0
	for acc in genomes:

		fileNumber+=1
		secondAminoAnalysis(genomes,acc,numGenomes,fileNumber,secondAmino,secondAminoTable4)


	# Close both of the output files
	secondAmino.close()
	secondAminoTable4.close()

#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
