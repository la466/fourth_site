#!/usr/bin/python

# Script number: 			4.5
# File: 					5 of 6
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Look at A usage over cds 5' region
# Output file(s):			_synonymous_sites.csv, _synonymous_sites_t4.csv, _nonsynonymous_sites.csv, _nonsynonymous_sites_t4.csv

import sys, os, csv, re

##########################
# VARIABLES
##########################

# Set the directory containing the sorted and parsed genomes
genomesDir = 'genome_extractions/cds/'
filteredGenesFile = "outputs/gene_filtering/_filteredGenes.csv"

testFiles = ["AP009152", "AP009153", "CP000817", "CP000925", "CP000264"]
testFile = ["AE005174"]

# Define the bases
bases = ["A", "C", "T", "G"]

codonNo = range(2,121)

output_dir = "outputs/ratio_testing/site_comparison/"

##########################
# FUNCTIONS
##########################

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

				# if acc in testFiles:
				genes.append(genus)
				genes.append(row[2:])
				genomes[acc] = genes


		numGenomes = len(genomes)

		return genomes, numGenomes

# Set up the new directories for the good files and the bad files
def setupDirectories(direc):

	# If the directory doesnt exist, make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)


# Define and open the output files
def newFile():

	# Output file with genes passing the filter
	output_fileFile = output_dir + "_synonymous_sites.csv"
	output_file = open(output_fileFile, "w")

	# Output file containing the genes failing the filter
	output_file_t4File = output_dir + "_synonymous_sites_t4.csv"
	output_file_t4 = open(output_file_t4File, "w")

	nonsynonymous_sites_path = output_dir + "_nonsynonymous_sites.csv"
	output_file_2 = open(nonsynonymous_sites_path, "w")

	nonsynonymous_sites_path_t4 = output_dir + "_nonsynonymous_sites_t4.csv"
	output_file_2_t4 = open(nonsynonymous_sites_path_t4, "w")

	return output_file, output_file_t4, output_file_2, output_file_2_t4


# Set up the file heading for the good genes
def fileHead(file):
	summaryHeadMisc = "Genus,Acc,trans_table,Num_Genes,GC3"

	for codon in codonNo:
		summaryHeadMisc += ",codon_%s" % (codon)

	summaryHeadMisc += "\n"

	file.write(summaryHeadMisc)




# Determine the proportion given a number and the total number
def proportion(number, total):
	prop = number / float(total)
	return prop



def output_fileAnalysis(genomes,acc,numGenomes,fileNumber,output_file,output_file_t4, output_file_2, output_file_2_t4):

	print '-' * 10
	print acc
	print "Genome %d of %d" % (fileNumber, numGenomes)

	genomeDir = genomesDir + acc + ".txt"

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

		fsite = {}
		# propSynsite = {}

		# for codon in codonNo:
		# 	fsite[codon] = {}
		# 	fsite[codon]["A"] = 0
		# 	fsite[codon]["T"] = 0
		# 	fsite[codon]["G"] = 0
		# 	fsite[codon]["C"] = 0

		# 	propSynsite[codon] = {}
		# 	for base in bases:
		# 		propSynsite[codon][base] = 0

		codonCount = 0
		codons = {}
		codons["A"] = 0
		codons["T"] = 0
		codons["G"] = 0
		codons['C'] = 0
		GC3 = 0

		gc_count = 0
		nt_count = 0


		synonymous_sites = {}
		for position in codonNo:
			synonymous_sites[position] = {}
			fsite[position] = {}
			for base in bases:
				synonymous_sites[position][base] = 0
				fsite[position][base] = 0

		gene_suitable_length = 0
		codons_suitable_length = 0


		for singleGene in genes:
			locus = re.findall('locus_tag=(.+?);', singleGene)[0]

			# Retrieve the gene sequence from the file
			seq = re.findall('\n(.+)', singleGene)[0]

			# Retrieve the gene sequence from the file
			getInfo["trans_table"] = re.findall('trans_table=(\d+?);', singleGene)[0]


			# If the locus is in the filtered loci that passed filtering
			if locus in loci:

				geneNumber += 1

				# Nt count (excluding stop)
				nt_count += len(seq) - 3

				# Get GC stats
				gc = re.subn('[GC]', '[GC]', seq[:-3])
				gc_count += gc[1]

				# Get GC details
				for i in range(0,len(seq)-3,3):
		 			codon = seq[i:i+3]

		 			# Get codon use (first site)
		 			if codon[0] in bases:
		 				codons[codon[0]] += 1
		 				codonCount+=1

		 			if codon[2] == "G" or codon[2] == "C":
		 				GC3 += 1

		 		# If the sequence is long enough
		 		if len(seq) > (max(codonNo)*3):

		 			gene_suitable_length += 1

		 			# For each codon, get the first and third codon positio
		 			for codon in codonNo:
		 				synonymous_site_seq_position = (codon*3) - 1
		 				synonymous_sites[codon][seq[synonymous_site_seq_position]] += 1

		 				fsite_seq_position = (codon*3)-3
		 				fsite[codon][seq[fsite_seq_position]] += 1


		# Calcalate the A proportions for each codon
		codon_prop = {}
		for codon in synonymous_sites:
			codon_prop[codon] = proportion(synonymous_sites[codon]["A"], gene_suitable_length)

		fsite_prop = {}
		for codon in fsite:
			fsite_prop[codon] = proportion(fsite[codon]["A"], gene_suitable_length)




		GC3Aprop = float(codons["A"])/codonCount
		GC3value = proportion(GC3, codonCount)

		gc_value = proportion(gc_count, nt_count)

		fileLine = "%s,%s,%s,%s,%s" % (genome,acc,getInfo["trans_table"],geneNumber,GC3value)
		fileLine2 = "%s,%s,%s,%s,%s" % (genome,acc,getInfo["trans_table"],geneNumber,gc_count)

		for codon in codon_prop:
			fileLine += ",%s" % codon_prop[codon]

		for codon in fsite_prop:
			fileLine2 += ",%s" % fsite_prop[codon]


		fileLine += "\n"
		fileLine2 += "\n"



		if getInfo["trans_table"] == "11":
			output_file.write(fileLine)
			output_file_2.write(fileLine2)
		else:
			output_file_t4.write(fileLine)
			output_file_2_t4.write(fileLine2)



#################################

def main():

	setupDirectories(output_dir)

	# Set the number of files to 0
	fileNumber = 0

	# Return in the information on the genomes
	genomes, numGenomes = getFilteredGenes()


	# Open the new files
	output_file,output_file_t4, output_file_2, output_file_2_t4 = newFile()

	# Write the heading for the new files

	file_head_1 = "Genus,Acc,trans_table,Num_Genes,GC3"
	file_head_2 = "Genus,Acc,trans_table,Num_Genes,gc"
	for codon in codonNo:
		file_head_1 += ",codon_%s" % (codon)
		file_head_2 += ",codon_%s" % (codon)

	file_head_1 += "\n"
	file_head_2 += "\n"

	output_file.write(file_head_1)
	output_file_t4.write(file_head_1)
	output_file_2.write(file_head_2)
	output_file_2_t4.write(file_head_2)

	fileNumber = 0
	for acc in genomes:

		#if acc in testFiles:
			fileNumber+=1
			output_fileAnalysis(genomes,acc,numGenomes,fileNumber,output_file,output_file_t4, output_file_2, output_file_2_t4)


	# Close both of the output files
	output_file.close()
	output_file_t4.close()
	output_file_2.close()
	output_file_2_t4.close()

#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
