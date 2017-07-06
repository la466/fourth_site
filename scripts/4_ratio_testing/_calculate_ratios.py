#!/usr/bin/python

# Script number: 			4.1
# File: 					1 of 6
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Calculated the ratios for each nucleotide at the given position
# Output file(s):			_site_SITE_ratios.csv, _startCodons.csv

import sys, os, csv, re, argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description = "Help for ratio calculations", formatter_class=RawTextHelpFormatter)
parser.add_argument("-site", help="Enter the site you wish to calculate ratios for.")
args = parser.parse_args()


if args.site:

	if float(args.site) > 15 or float(args.site) < 0:
		print "\nWARNING: You have chosen an invalid site. Choose a site between 4 and 15, or 0 for all sites between 5 and 15 .\n"
		sys.exit()
	else:
		if int(args.site) == 0:
			print "\nCalculating ratios for site all sites\n"

else:
	sys.exit()

##########################
# VARIABLES
##########################


# Set the directory containing the sorted and parsed genomes
genomesDir = 'genome_extractions/cds/'

output_directory = "outputs/ratio_testing/"

# site_directory = output_directory + "site_%s/" % args.site

# testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
# testFile = ["AE005174"]

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
def newFile(site):


	site_directory = output_directory + "site_%s/" % site
	setupDirectories(site_directory)


	# Output file with genes passing the filter
	fourthSiteFile = site_directory + "_site_%d_ratios.csv" % site
	fourthSite = open(fourthSiteFile, "w")

	# Output file containing the genes failing the filter
	controlFile = site_directory + "_site_%d_ratios_t4.csv" % site
	control = open(controlFile, "w")

	absApath = site_directory + "_site_%d_abs_A_count.csv" % site
	absAfile = open(absApath, "w")

	absApathT4 = site_directory + "_site_%d_abs_A_count_t4.csv" % site
	absAfileT4 = open(absApathT4, "w")

	startCodonPath = site_directory + "_site_%d_start_codon_A_ratios.csv" % site
	startCodonFile = open(startCodonPath, "w")

	startCodonPathT4 = site_directory + "_site_%d_start_codon_A_ratios_t4.csv" % site
	startCodonFileT4 = open(startCodonPathT4, "w")

	return fourthSite, control, absAfile, absAfileT4, startCodonFile, startCodonFileT4


# Set up the file heading for the good genes
def fileHead(file, site):



	head_line = "Genus,Acc,trans_table,Num_Genes,GC,GC3,"

	for base in bases:
		summaryHead = "Prop_%s%d,Obs_%s%d,%s_Codons,Prop_%s_Codons,%s_Ratio," % (base, site, base, site, base, base, base)
		head_line += summaryHead

	head_line += "Codon_Count"
	head_line += "\n"
	file.write(head_line)



def fileHead2(file):

	headLine = "genus,acc,trans_table,num_genes,gc"
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



def fourthSiteAnalysis(genomes,acc,numGenomes,fileNumber,fourthSite,control, absAfile, absAfileT4, startCodonFile, startCodonFileT4, site):

	print '-' * 10
	print acc
	print "Genome %d of %d" % (fileNumber, numGenomes)

	genomeDir = 'genome_extractions/cds/' + acc + ".txt"

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
			phylum = re.findall('phylum=(.+?);', singleGene)[0]



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
		 		pos4 = seq[site-1]
		 		if pos4 in bases:
		 			baseFour[pos4] += 1

		 		# Get the prortion of codons using nucletoides
		 		for i in range(0,len(seq)-3,3):
		 			codon = seq[i:i+3]


		 			if site % 3 == 1:
		 				codonPos = 0
		 			elif site % 3 == 2:
		 				codonPos = 1
		 			elif site % 3 == 0:
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
		summaryA = "%s,%s,%s,%s,%s," % (propA,fourthA, codonA, propCodonA, propA/propCodonA)
		summaryT = "%s,%s,%s,%s,%s," % (propT,fourthT, codonT, propCodonT, propT/propCodonT)
		summaryG = "%s,%s,%s,%s,%s," % (propG,fourthG, codonG, propCodonG, propG/propCodonG)
		summaryC = "%s,%s,%s,%s,%s," % (propC,fourthC, codonC, propCodonC, propC/propCodonC)
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
			fourthSite.write(summaryMisc + summaryA + summaryT + summaryG + summaryC + summaryEnd)
			absAfile.write(absALine)
			startCodonFile.write(startCodonLine)
		else:
			control.write(summaryMisc + summaryA + summaryC + summaryT + summaryG + summaryEnd)
			absAfileT4.write(absALine)
			startCodonFileT4.write(startCodonLine)



#################################

def main():

	setupDirectories(output_directory)

	if int(args.site) == 0:
		sites = [4,5,6,7,8,9,10,11,12,13,14,15]
	elif int(args.site) != 0:
		sites = [args.site]




	for site in sites:


		site = int(site)



		print "\nCalculating ratios for site %s\n" % site

		# Set the number of files to 0
		fileNumber = 0

		# Return in the information on the genomes
		genomes, numGenomes = getFilteredGenes()


		# Open the new files
		fourthSite,control, absAfile, absAfileT4, startCodonFile, startCodonFileT4 = newFile(int(site))

		# Write the heading for the new files
		fileHead(fourthSite, int(site))
		fileHead(control, int(site))

		fileHead2(absAfile)
		fileHead2(absAfileT4)

		# Set up the start codon file
		headLine = "acc"
		for base in ["a", "c", "g", "t"]:
			headLine += ",%stg_a%d_count,%stg_count,%stg_a%d_prop,%s_codon_prop,%stg_a%d_ratio" % (base, site, base, base, site,base,base, site)
		headLine += "\n"
		startCodonFile.write(headLine)
		startCodonFileT4.write(headLine)

		fileNumber = 0
		for acc in genomes:
			# if acc in testFile:
				fileNumber+=1
				fourthSiteAnalysis(genomes,acc,numGenomes,fileNumber,fourthSite,control,absAfile,absAfileT4, startCodonFile, startCodonFileT4, int(site))


		# Close the output files
		fourthSite.close()
		control.close()
		absAfile.close()
		absAfileT4.close()

#################################


if __name__ == "__main__":
	main()
