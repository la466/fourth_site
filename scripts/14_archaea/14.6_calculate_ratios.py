#!/usr/bin/python

# Script number: 			14.6
# File: 					6 of 7
# Prerequisite script(s):	14.1, 14.2, 14.3, 14.4, 14.5
# Prerequisite file(s):		_filteredGenes_archaea.csv
# Description: 				Calculated the ratios for each nucleotide at the given position
# Output file(s):			_site_SITE_ratios.csv, _startCodons.csv, table4genomes.txt

import sys, os, csv, re, argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description = "Help for ratio calculations", formatter_class=RawTextHelpFormatter)
parser.add_argument("site", help="Enter the site you wish to calculate ratios for.")
args = parser.parse_args()

args.site = int(args.site)

if args.site:

	if float(args.site) > 12 or float(args.site) < 0:
		print "\nWARNING: You have chosen a site greater than 12. Please re-run script 14.6 (14.6_calculate_ratios.py) using another site.\n"
		sys.exit()
	else:
		print "\nCalculating ratios for site %d\n" % args.site
else:
	sys.exit()

##########################
# VARIABLES
##########################


# Set the directory containing the sorted and parsed genomes
genomesDir = 'archaea_gneomes/cds/'

output_directory = "outputs/archaea/ratio_testing/"

site_directory = output_directory + "site_%s/" % args.site

testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE005174"]

# Define the bases
bases = ["A", "T", "G", "C"]



##########################
# FUNCTIONS
##########################

# Set up the new directories for the good files and the bad files
def setupDirectories(direc):

	# If the directory doesnt exist, make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)


# Get all the good genes
def getFilteredGenes():

	goodGenes = {}

	filteredGenesFile = "outputs/archaea/gene_filtering/_filteredGenes_archaea.csv"
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
	fourthSiteFile = site_directory + "_site_%d_ratios_archaea.csv" % args.site
	fourthSite = open(fourthSiteFile, "w")

	# Output file containing the genes failing the filter
	controlFile = site_directory + "_site_%d_ratios_t4_archaea.csv" % args.site
	control = open(controlFile, "w")


	return fourthSite, control


# Set up the file heading for the good genes
def fileHead(file):

	head_line = "Genus,Acc,trans_table,Num_Genes,GC,GC3,"

	for base in bases:
		summaryHead = "Prop_%s%d,Obs_%s%d,%s_Codons,Prop_%s_Codons,%s_Ratio," % (base, args.site, base, args.site, base, base, base)
		head_line += summaryHead

	head_line += "Codon_Count"
	head_line += "\n"
	file.write(head_line)



def fileHead2(file):

	headLine = "genus,phylum,acc,trans_table,num_genes,gc"
	headLine += ",a_second,prop_a,a_codons,a_codon_prop"
	headLine += ",c_second,prop_c,c_codons,c_codon_prop"
	headLine += ",t_second,prop_t,t_codons,t_codon_prop"
	headLine += ",g_second,prop_g,g_codons,g_codon_prop"
	headLine += "\n"

	file.write(headLine)


# Determine the proportion given a number and the total number
def proportion(number, total):

	if total != 0:
		prop = number / float(total)
	else:
		prop = 0
	return prop



def fourthSiteAnalysis(genomes,acc,numGenomes,fileNumber,fourthSite,control):

	print '-' * 10
	print acc
	print "Genome %d of %d" % (fileNumber, numGenomes)

	genomeDir = 'archaea_genomes/cds/' + acc + ".txt"

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
		getInfo["phylum"] = ""

		baseFour = {}
		baseFour["A"] = 0
		baseFour["T"] = 0
		baseFour["G"] = 0
		baseFour["C"] = 0

		codonCount = 0
		codons = {}
		codons["A"] = 0
		codons["T"] = 0
		codons["G"] = 0
		codons['C'] = 0

		GC3 = 0
		GC = 0
		nucleotides = 0

		# Set up the start codon count
		startCodonCount = {}
		for base in bases:
			startCodonCount[base] = []

		# For each of the cds
		for singleGene in genes:
			locus = re.findall('locus_tag=(.+?);', singleGene)[0]


			# Retrieve the gene sequence from the file
			seq = re.findall('\n(.+)', singleGene)[0]

			# Get the GC content of the cds
			gc = re.subn('[GC]', "0", seq[:-3])
			GC += gc[1]
			nucleotides += len(seq[:-3])

			# Retrieve the gene sequence from the file
			getInfo["trans_table"] = re.findall('trans_table=(\d+?);', singleGene)[0]


			# If the locus is in the filtered loci that passed filtering
			if locus in loci:

				geneNumber += 1


				# Get the nucleotide at the queried position (pos4 was original)
		 		pos4 = seq[args.site-1]
		 		if pos4 in bases:
		 			baseFour[pos4] += 1

		 		# Get the prortion of codons using nucletoides
		 		for i in range(0,len(seq)-3,3):
		 			codon = seq[i:i+3]


		 			if args.site % 3 == 1:
		 				codonPos = 0
		 			elif args.site % 3 == 2:
		 				codonPos = 1
		 			elif args.site % 3 == 0:
		 				codonPos = 2


		 			if codon[codonPos] in bases:
		 				codons[codon[codonPos]] += 1
		 				codonCount+=1


		 			# Get GC3 content
		 			if codon[2] == "G" or codon[2] == "C":
		 				GC3 += 1

				if pos4 == "A":
					startCodonCount[seq[0]].append(1)
				else:
					startCodonCount[seq[0]].append(0)


		# Genome calculations
		GCcontent = proportion(GC, nucleotides)

		fourthA = baseFour["A"]
		fourthT = baseFour["T"]
		fourthG = baseFour["G"]
		fourthC = baseFour["C"]

		propA = proportion(fourthA,geneNumber)
		propT = proportion(fourthT,geneNumber)
		propG = proportion(fourthG,geneNumber)
		propC = proportion(fourthC,geneNumber)




		codonA = codons["A"]
		codonT = codons["T"]
		codonG = codons["G"]
		codonC = codons["C"]


		propCodonA = proportion(codonA,codonCount)
		propCodonT = proportion(codonT,codonCount)
		propCodonG = proportion(codonG,codonCount)
		propCodonC = proportion(codonC,codonCount)




		# File lines
		summaryMisc = "%s,%s,%s,%s,%s,%s," % (genome,acc,getInfo["trans_table"],geneNumber,GCcontent,proportion(GC3,codonCount))
		summaryA = "%s,%s,%s,%s,%s," % (propA,fourthA, codonA, propCodonA, proportion(propA, propCodonA))
		summaryT = "%s,%s,%s,%s,%s," % (propT,fourthT, codonT, propCodonT, proportion(propT, propCodonT))
		summaryG = "%s,%s,%s,%s,%s," % (propG,fourthG, codonG, propCodonG, proportion(propG, propCodonG))
		summaryC = "%s,%s,%s,%s,%s," % (propC,fourthC, codonC, propCodonC, proportion(propC, propCodonC))
		summaryEnd = "%s\n" % codonCount

		absALine = "%s,%s,%s,%s,%s," % (genome,acc,getInfo["trans_table"],geneNumber,GCcontent)
		absALine += "%s,%s,%s,%s," % (fourthA, propA, codonA, propCodonA)
		absALine += "%s,%s,%s,%s," % (fourthC, propC, codonC, propCodonC)
		absALine += "%s,%s,%s,%s," % (fourthT, propT, codonT, propCodonT)
		absALine += "%s,%s,%s,%s" % (fourthG, propG, codonG, propCodonG)
		absALine += "\n"


		startCodonLine = "%s" % acc
		for base in startCodonCount:
			propA = proportion(sum(startCodonCount[base]),len(startCodonCount[base]))
			AcodonProp = proportion(codons[base], codonCount)
			startCodonLine += ",%s,%s,%s,%s,%s" % (sum(startCodonCount[base]),len(startCodonCount[base]),propA, AcodonProp, proportion(propA, AcodonProp))
		startCodonLine += "\n"



		# Write to the files
		if getInfo["trans_table"] == "11":
			if geneNumber != 0:

				fourthSite.write(summaryMisc + summaryA + summaryT + summaryG + summaryC + summaryEnd)

		else:
			if geneNumber != 0:
				control.write(summaryMisc + summaryA + summaryC + summaryT + summaryG + summaryEnd)




#################################

def main():

	setupDirectories(output_directory)
	setupDirectories(site_directory)

	# Set the number of files to 0
	fileNumber = 0

	# Return in the information on the genomes
	genomes, numGenomes = getFilteredGenes()


	# Open the new files
	fourthSite,control = newFile()

	# Write the heading for the new files
	fileHead(fourthSite)
	fileHead(control)



	fileNumber = 0
	for acc in genomes:
		fileNumber+=1
		if fileNumber:
			fourthSiteAnalysis(genomes,acc,numGenomes,fileNumber,fourthSite,control)


	# Close the output files
	fourthSite.close()
	control.close()


#################################


if __name__ == "__main__":
	main()
