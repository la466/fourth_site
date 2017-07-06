#!/usr/bin/python

# Script number: 			4.3
# File: 					3 of 6
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Calculate ratios having removed 4 site overlapping genes
# Output file(s):			_site_4_no_overlap.csv, _site_4_no_overlap_stats.txt

import sys, os, csv, re


##########################
# VARIABLES
##########################


# Set the directory containing the sorted and parsed genomes
genomesDir = 'bacterial_genomes/codons/genomes/'

testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE005174"]

# Define the bases
bases = ["A", "T", "G", "C"]


overlap_stats_path = "outputs/ratio_testing/_sie_4_no_overlap_stats.txt"

##########################
# FUNCTIONS
##########################


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

	outputDir = "outputs/ratio_testing/site_4/"

	if not os.path.exists(outputDir):
		os.makedirs(outputDir)


	# Output file with genes passing the filter
	fourthSiteFile = outputDir + "_site_4_no_overlap.csv"
	fourthSite = open(fourthSiteFile, "w")



	return fourthSite


# Set up the file heading for the good genes
def fileHead(file):
	summaryHeadMisc = "Genus,Acc,trans_table,Num_Genes,GC3,"
	summaryHeadA = "Prop_A4,Obs_A4,A_Codons,Prop_A_Codons,A_Ratio,"
	summaryHeadT = "Prop_T4,Obs_T4,T_Codons,Prop_T_Codons,T_Ratio,"
	summaryHeadG = "Prop_G4,Obs_G4,G_Codons,Prop_G_Codons,G_Ratio,"
	summaryHeadC = "Prop_C4,Obs_C4,C_Codons,Prop_C_Codons,C_Ratio,"
	summaryHeadP = "Codon_Count"
	summaryHeadBreak = "\n"
	file.write(summaryHeadMisc+summaryHeadA+summaryHeadC+summaryHeadT+summaryHeadG+summaryHeadP+summaryHeadBreak)


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
	prop = number / float(total)
	return prop



def fourthSiteAnalysis(genomes,acc,numGenomes,fileNumber,fourthSite, total_genes, total_4_overlap):

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

		seq_4_a = 0


		GC3 = 0
		GC = 0
		nucleotides = 0

		overlaps = {}
		overlaps[0] = {}
		overlaps[1] = {}

		overlapLeading = 0
		overlapLagging = 0

		nonOverlaps = []
		overlapLoci = []

		total_cds_count = 0


		for singleGene in genes:

			locus = re.findall('locus_tag=(.+?);', singleGene)[0]

			# If the gnee has passed filtering
			if locus in loci:

				total_cds_count += 1

				position = re.findall('loc=(.+?);', singleGene)[0]
				strand = int(re.findall('complement=(\d+?);', singleGene)[0])
				seq = re.findall('\n(.+)', singleGene)[0]

				# gc = re.subn('[GC]', "", seq[:-3])
				# GC += gc[1]
				# nucleotides += len(seq[:-3])

				if seq[3] == "A":
					seq_4_a += 1

				# for i in range(0,len(seq)-3,3):
			 # 		codon = seq[i:i+3]
			 # 		if codon[2] == "G" or codon[2] == "C":
			 # 			GC3 += 1

				cdsStart = re.findall('^(\d+)(?=\.\.)', position)[0]
				cdsEnd = re.findall('(?<=\.\.)(\d+)$', position)[0]



				# If it he leading strand, add to array of leading
				if strand == 0:
					overlapLeading += 1
					overlaps[0][overlapLeading] = [int(cdsStart), int(cdsEnd), locus, seq]

				# If the lagging strand, add to the array of lagging
				elif strand == 1:
					overlapLagging += 1
					overlaps[1][overlapLagging] = [int(cdsStart), int(cdsEnd), locus, seq]


		# For each of the leading and lagging strands
		for strand in overlaps:

				# For each of the cds on those strands
				for cdsID in overlaps[strand]:

					# If the cds is on the leading strand
					if strand == 0:

						# If is not the first cds (no overlap)
						if cdsID != 1:

							focusStart = overlaps[strand][cdsID][0]
							prevEnd = overlaps[strand][cdsID-1][1]

							if prevEnd - focusStart == 3:
								overlapLoci.append(overlaps[strand][cdsID][2])


					elif strand == 1:

						# If the cds is not the first one the lagging strand
						if cdsID != len(overlaps[strand]):

							focusEnd = overlaps[strand][cdsID][1]
							focusStart = overlaps[strand][cdsID][0]
							prevEnd = overlaps[strand][cdsID+1][1]
							prevStart = overlaps[strand][cdsID+1][0]

							if focusEnd - prevStart == 3:
								overlapLoci.append(overlaps[strand][cdsID][2])




		# For each gene
		for singleGene in genes:
			locus = re.findall('locus_tag=(.+?);', singleGene)[0]

			# Retrieve the gene sequence from the file
			seq = re.findall('\n(.+)', singleGene)[0]

			# If the gene passed filtering
			if locus in loci and locus not in overlapLoci:

				# Retrieve the gene sequence from the file
				getInfo["trans_table"] = re.findall('trans_table=(\d+?);', singleGene)[0]


				for i in range(0,len(seq)-3,3):
			 		codon = seq[i:i+3]


				gc = re.subn('[GC]', "0", seq[:-3])
				GC += gc[1]
				nucleotides += len(seq[:-3])

				geneNumber += 1
		 		pos4 = seq[3]
		 		if pos4 in bases:
		 			baseFour[pos4] += 1

		 		for i in range(0,len(seq)-3,3):
		 			codon = seq[i:i+3]



		 			if codon[0] in bases:
		 				codons[codon[0]] += 1
		 				codonCount+=1



		 			if codon[2] == "G" or codon[2] == "C":
		 				GC3 += 1



		total_genes.append(total_cds_count)
		total_4_overlap.append(len(overlapLoci))

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

		# codonCount = numpy.sum([codonA,codonT,codonG,codonC])
		propCodonA = proportion(codonA,codonCount)
		propCodonT = proportion(codonT,codonCount)
		propCodonG = proportion(codonG,codonCount)
		propCodonC = proportion(codonC,codonCount)

		summaryMisc = "%s,%s,%s,%s,%s," % (genome,acc,getInfo["trans_table"],geneNumber,proportion(GC3,codonCount))
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

		if getInfo["trans_table"] == "11":
			fourthSite.write(summaryMisc + summaryA + summaryC + summaryT + summaryG + summaryEnd)





#################################

def main():

	# Set the number of files to 0
	fileNumber = 0

	# Return in the information on the genomes
	genomes, numGenomes = getFilteredGenes()


	# Open the new files
	fourthSite = newFile()

	# Write the heading for the new files
	fileHead(fourthSite)


	total_4_overlap = []
	total_genes = []



	fileNumber = 0
	for acc in genomes:
		#if acc in testFile:
			fileNumber+=1
			fourthSiteAnalysis(genomes,acc,numGenomes,fileNumber,fourthSite, total_genes, total_4_overlap)

	overlap_percentage = sum(total_4_overlap)/float(sum(total_genes))
	overlap_stats = open(overlap_stats_path, "w")
	overlap_line = "Total 4 site overlaps: %d\n" % sum(total_4_overlap)
	overlap_line += "Total cds: %d\n" % sum(total_genes)
	overlap_line += "4 site overlap prop: %s" %	overlap_percentage
	overlap_stats.write(overlap_line)

	# Close both of the output files
	fourthSite.close()
	overlap_stats.close()


#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
