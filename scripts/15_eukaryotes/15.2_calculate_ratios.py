#!/usr/bin/python

# Script number: 			15.2
# File: 					2 of 2
# Prerequisite script(s):	15.1
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Calculated the ratios for each nucleotide at the given position
# Output file(s):			_site_SITE_ratios.csv, _startCodons.csv, table4genomes.txt

import sys, os, csv, re, argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO

parser = argparse.ArgumentParser(description = "Help for ratio calculations", formatter_class=RawTextHelpFormatter)
parser.add_argument("-eukaryote_site", help="Enter the site you wish to calculate ratios for.")
parser.add_argument("-genomes", help="Enter the eukaroytic genome you wish to calculate ratios for. Enter eukaryotic filename or 'all_genomes'.")
args = parser.parse_args()


if args.genomes:

	if args.genomes == "all_genomes":
		genome_select = "all_genomes"
	else:
		genome_select = args.genomes
else:	
	print "\nPlease select a genome to calculate ratios for all pass 'all_genomes'.\n"
	sys.exit()


if args.eukaryote_site:

	if type(args.eukaryote_site) == "int":
		if float(args.eukaryote_site) > 12 or float(args.eukaryote_site) < 0:
			print "\nWARNING: You have chosen a site greater than 12. Please re-run script 15.2 (15.2_calculate_ratios.py) using another site.\n"
			sys.exit()
		else:
			print "\nCalculating ratios for site %d\n" % args.eukaryote_site
else:
	sys.exit()



##########################
# VARIABLES
##########################


sites = [4,5,6,7,8,9,10,11,12]

# Set the directory containing the sorted and parsed genomes
genomesDir = 'eukaryote_genomes/'

output_directory = "outputs/eukaryotes/"

#site_directory = output_directory + "site_%s/" % args.eukaryote_site

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

	filteredGenesFile = "outputs/eukaryotes/gene_filtering/_filteredGenes_eukaryotes.csv"
	with open(filteredGenesFile, 'U') as myfile:
		read = csv.reader(myfile)
		
		genomes = {}

		rowNum = 0
		for row in read:		
			rowNum += 1
			if rowNum > 1:
				

				genus = row[0]
				acc = row[1]
				
				#print genus

				genes = []

				# if acc in testFile: 
				genes.append(genus)
				genes.append(row[2:])
				genomes[genus] = genes

			
		numGenomes = len(genomes)


		return genomes, numGenomes




# Set up the file heading for the good genes
def fileHead(file, site):	

	head_line = "Genus,Num_Genes,GC,GC3,"

	for base in bases:
		summaryHead = "Prop_%s%d,Obs_%s%d,%s_Codons,Prop_%s_Codons,%s_Ratio," % (base, site, base, site, base, base, base)
		head_line += summaryHead
	
	head_line += "Codon_Count"
	head_line += "\n"
	file.write(head_line)





# Determine the proportion given a number and the total number
def proportion(number, total):

	if total != 0:
		prop = number / float(total)
	else: 
		prop = 0
	return prop



def fourthSiteAnalysis(genomes,genome,acc,numGenomes,fileNumber,fourthSite, geneNumber, getInfo, baseFour, codonCount, codons,GC3, GC, nucleotides, startCodonCount, site, filtered_genes):

	print '-' * 10
	print acc
	print "Genome %d of %d" % (fileNumber, numGenomes)
	print "Calculating ratios for site %s" % site

	genomeDir = genomesDir + genome
	print "\n%s\n" % genomeDir
	
	genome = genomes[acc][0]
	loci = filtered_genes

	cds_count = 0
	codon_count = 0
	GC = 0
	GC3 = 0
	nucleotide_count = 0


	# with open(genomeDir) as myfile:
	# 	line = myfile.read()
	# 	genes = line.split("\n\n")
	# 	genes.pop()

	# 	# Set the gene count for the genome to 0


	# 	# For each of the cds
	# 	for singleGene in genes:

	handle = open(genomeDir, "rU")
	for record in SeqIO.parse(handle, "fasta"):

		#locus = re.findall('locus_tag=(.+?);', singleGene)[0]
		locus = record.id
                    
		# Retrieve the gene sequence from the file
		#seq = re.findall('\n(.+)', singleGene)[0]
		seq = str(record.seq)

		# Get the GC content of the cds
		gc = re.subn('[GC]', "0", seq[:-3])
		GC += gc[1]
		nucleotide_count += len(seq[:-3])

		# Retrieve the gene sequence from the file
		
		

		# If the locus is in the filtered loci that passed filtering
		if locus in loci:

			cds_count += 1

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
	 				codon_count+=1


	 			# Get GC3 content
	 			if codon[2] == "G" or codon[2] == "C":
	 				GC3 += 1
			
			if pos4 == "A":
				startCodonCount[seq[0]].append(1)
			else:
				startCodonCount[seq[0]].append(0)


	return cds_count, GC, GC3, codon_count, nucleotide_count


def get_genomes(direc):

	genome_list = []

	# For each file in the directory
	for eachfile in os.listdir(direc):

		if eachfile.endswith('.fa'):

			genome_list.append(eachfile)

	return genome_list


#################################

def main():

	# Get an array of the genomes to be tested
	if genome_select == "all_genomes":
		genome_list = get_genomes(genomesDir)
	else:
		
		genome_list = []
		genome_list.append(genome_select)





	# Setup the directories to contain the outputs
	for genome in genome_list:


		genus_string = re.findall("^(.*?)\..*", genome)[0] 
		genus_string = genus_string.lower()

		genome_dir_path = output_directory + genus_string + "/"
		setupDirectories(genome_dir_path)

		if args.eukaryote_site == "all_sites":

			for site in sites:
				site_path = "%ssite_%s/" % (genome_dir_path, site)
				setupDirectories(site_path)

		else:
			site = str(args.eukaryote_site)
			site_path = genome_dir_path + "site_" + site + "/"
			setupDirectories(site_path)




	# Get an array of sites to query for the genome
	site_query = []
	if args.eukaryote_site == "all_sites":
		site_query = sites
	else:
		site_query.append(int(args.eukaryote_site))


	# Return in the information on the genomes
	genomes, numGenomes = getFilteredGenes()

	numGenomes = len(genome_list)

	
	fileNumber = 0
	for genome in genome_list:


		fileNumber += 1
		

		genus = re.findall("^(.*?)\..*", genome)[0] 
		genus = genus.lower()
		

		filtered_genes = genomes[genus][1]
		


		for site in site_query:


	
			# Open the new files
			genome_site_path = "outputs/eukaryotes/%s/site_%s/_site_%s_ratios.csv" % (genus, site, site)
			fourthSite = open(genome_site_path, "w")

			# Write the heading for the new files
			fileHead(fourthSite, site)




			geneNumber = 0
			getInfo = {}

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



			
			acc = genus
			cds_count, gc_count, gc3_count, codon_count, nucleotide_count = fourthSiteAnalysis(genomes,genome,acc,numGenomes,fileNumber,fourthSite, geneNumber, getInfo, baseFour, codonCount, codons,GC3, GC, nucleotides, startCodonCount, site, filtered_genes)
			

			geneNumber += cds_count
			GC += gc_count
			GC3 += gc3_count
			codonCount += codon_count
			nucleotides += nucleotide_count

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
			summaryMisc = "%s,%s,%s,%s," % (genus,geneNumber,GCcontent,proportion(GC3,codonCount))
			summaryA = "%s,%s,%s,%s,%s," % (propA,fourthA, codonA, propCodonA, proportion(propA, propCodonA))
			summaryT = "%s,%s,%s,%s,%s," % (propT,fourthT, codonT, propCodonT, proportion(propT, propCodonT))
			summaryG = "%s,%s,%s,%s,%s," % (propG,fourthG, codonG, propCodonG, proportion(propG, propCodonG))
			summaryC = "%s,%s,%s,%s,%s," % (propC,fourthC, codonC, propCodonC, proportion(propC, propCodonC))
			summaryEnd = "%s\n" % codonCount




			# Write to the files
			fourthSite.write(summaryMisc + summaryA + summaryT + summaryG + summaryC + summaryEnd)




			# Close the output files
			fourthSite.close()

#################################


if __name__ == "__main__":
	main()
