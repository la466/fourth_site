#!/usr/bin/python

# Script number: 			8.5
# File: 					5 of 5
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 4.1, 8.1, 8.2, 8.3
# Prerequisite file(s):		ACC/ACC.out, _site_4_ratios.csv
# Description: 				Analyse the CAI values from CodonW run.
# Output file(s):			start_codon_CAI_means.csv


##########################
# VARIABLES
##########################




import sys, os, re, shutil, Bio, csv
import numpy as np
from sys import argv


genomes_dir = "genome_extractions/cds/"
expression_dir = "outputs/expression/"
expression_genomes_dir = expression_dir + "genomes/"

testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE005174"]




##########################
# FUNCTIONS
##########################



# Create a list of all the good genomes
def getGenomes(direc):

	files = []

	# For each file in the downloaded genomes
	for eachfile in os.listdir(direc):

		if eachfile != '.DS_Store':
			# if eachfile in testFile:
				files.append(eachfile)

	return files, len(files)

# Split the text file containing the CDS for the genome
def splitFile(file, split):
	with open(file) as myfile:
		line = myfile.read()
		genes = line.split(split)

		return genes


def get_t4_genomes(file_path):

	t4_genomes = []

	with open(file_path, 'rU') as myfile:

		lines = myfile.readlines()

		for line in lines:
			line = line.strip('\n')
			splits = line.split(',')
			t4_genomes.append(splits[1])

	return(t4_genomes)



def proportion(number, total):

	if total != 0:
		pc = number / float(total)
	else:
		pc = 0
	return pc






def getCAIs(acc, CAIvals):



	CAIvals[acc] = {}


	# # Split the file containing the CAI values and
	# # remove the first and last entry (blank)
	CAIfilePath = expression_genomes_dir + acc + "/" + acc + ".out"
	CAIs = splitFile(CAIfilePath, "\n")
	CAIs.pop()

	CAIline = 0

	for cds in CAIs:

		CAIline +=1

		# For each of the CAIs
		if CAIline > 1:

			# Split the line
			cdsLine = cds.split("\t")

			# Retrieve the gene name and the CAI value
			geneName = cdsLine[0]
			CAIval = cdsLine[8]
			CAIval = float(CAIval)

			CAIvals[acc][geneName] = CAIval


def calcMeanCAIs(HE_genes, CAIvals, acc):

	meanCAIs = {}

	for acc in CAIvals:

		CAItotal = 0
		CAInumber = 0

		# Calculate the mean CAI
		for gene in CAIvals[acc]:

			#Ensure the gene isnt one of the refernece genes
			if gene not in HE_genes[acc]:

				CAItotal += CAIvals[acc][gene]

				CAInumber +=1

		if CAItotal == 0:
			meanCAI = 0
		else:
			meanCAI = CAItotal / float(CAInumber)

		meanCAIs[acc] = meanCAI

	return meanCAIs

def get_genome(fourthPath):


	GC3s = {}
	Aratios = {}

	with open(fourthPath, "U") as myfile:

		lineNum = 0
		for line in myfile:

			lineNum +=1

			if lineNum >1:
				splits = line.split(",")
				acc = splits[1]
				GC3 = splits[4]
				Aratio = splits[9]

				GC3s[acc] = GC3
				Aratios[acc] = Aratio

	return GC3s, Aratios




def write_mean_CAIGC(CAIGCpath, meanCAIs, GC3s, Aratios):


	CAIGC = open(CAIGCpath, "w")
	CAIGC.write("acc,gc3,meanCAI,a_ratios\n")

	for acc in meanCAIs:

		if acc in GC3s:

			if meanCAIs[acc] != 0:
				fileLine = "%s,%s,%s,%s\n" % (acc,GC3s[acc], meanCAIs[acc], Aratios[acc])
				CAIGC.write(fileLine)



def get_HE_genes(CAIvals):


	HEgenes = {}


	# Get a list of the highly expressed genes used
	for acc in CAIvals:

		HEgenes[acc] = []


		# Get the genes used as highly expressed genes for each genome
		hePath = expression_genomes_dir + acc + "/" + acc + "_he.txt"

		with open(hePath, "U") as myfile:

			line = myfile.read()
			genes = line.split("\n")


			for singleHEgene in genes:

				gene = re.findall('^>(.+)(?=\\t\\t\\t)', singleHEgene)



				if len(gene) != 0:
					HEgenes[acc].append(gene[0])

	return HEgenes



def get_geneSeq(acc, geneSeq):


	geneSeq[acc] = {}

	genomeFile = genomes_dir + acc + ".txt"

	with open(genomeFile, "U") as myfile:

		line = myfile.read()
		genes = line.split("\n\n")
		genes.pop()

		unknownCount = 0

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


			# print gene, len(seq)

			geneSeq[acc][gene] = seq

			# print acc, gene, geneSeq[acc][gene]

	return geneSeq


def start_codon_CAI(HE_genes, CAIvals, geneSeq):


	start_codons = {}

	# For each of the genomes with CAI values
	for acc in CAIvals:

		cds_count = 0
		start_codons[acc] = {}

		for start in ['ATG', 'GTG', 'TTG']:
			start_codons[acc][start] = []

		# For each of the gene loci in the genome
		for locus in CAIvals[acc]:

			# Ensure the gnee is not one of the reference genes
			if locus not in HE_genes[acc]:


				if CAIvals[acc][locus] != 0:

					# If that locus is also in the gene sequences
					if locus in geneSeq[acc]:

						CAIval = CAIvals[acc][locus]
						seq = geneSeq[acc][locus]

						cds_count += 1
						if seq[:3] in ['ATG', 'GTG', 'TTG']:
							start_codons[acc][seq[:3]].append(CAIval)



	return(start_codons)






#################################

def main():


	fileNumber = 0

	# Get a list of the genomes containing CAI information
	genomes, numGenomes = getGenomes(expression_genomes_dir)

	t4_genomes = get_t4_genomes('outputs/gene_filtering/table4genomes.txt')



	# Get the GC3 values for each genome
	print "\nGetting the GC3 values for each genome\n"
	# Get the seq for each gene
	print "\nGetting the gene sequence for each gene\n"



	fourth_site_path = "outputs/ratio_testing/site_4/_site_4_ratios.csv"
	GC3s, Aratios = get_genome(fourth_site_path)
	CAIvals = {}
	geneSeq = {}


	for acc in genomes:
		# if acc in testFile:
		fileNumber += 1

		if fileNumber:

			print "-" * 10
			print acc


			# Get the CAI values for each genome
			getCAIs(acc, CAIvals)
			geneSeq = get_geneSeq(acc, geneSeq)






	fileNumber = 0

	# Get a list of the highly expressed ribosomal genes used as
	# a reference set
	HE_genes = get_HE_genes(CAIvals)

	# Calculate the mean CAI for the genome
	# meanCAIs = calcMeanCAIs(HE_genes, CAIvals, acc)

	start_codons = start_codon_CAI(HE_genes, CAIvals, geneSeq)


	output_file = open('outputs/expression/start_codon_CAI_means.csv', 'w')
	output_file.write('acc,atg_cai,gtg_cai,ttg_cai\n')



	for acc in start_codons:

		if acc not in t4_genomes:

			output_line = '%s' % acc

			highest = ''
			highest_start = 0

			for start in sorted(start_codons[acc]):
				mean_cai = np.mean(start_codons[acc][start])
				output_line += ',%s' % (mean_cai)





			output_line += '\n'
			output_file.write(output_line)

	output_file.write(output_line)

	output_file.close()



#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
