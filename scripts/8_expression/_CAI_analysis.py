#!/usr/bin/python

# Script number: 			8.4
# File: 					4 of 5
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 4.1, 8.1, 8.2, 8.3
# Prerequisite file(s):		ACC/ACC.out, _site_4_ratios.csv
# Description: 				Analyse the CAI values from CodonW run.
# Output file(s):			_ecoli_codonw_compare.csv, _CAI_a_ratios.csv, _CAI_a_notA.csv, _CAI_a_nota_restrictedGC3.csv, _high_low_CAI_compare.csv, high_low_CAI_compare_restrictedGC3.csv, _lowMeanCAIcompareFourthA.csv, _highMeanCAIcompareFourthA.csv


##########################
# VARIABLES
##########################


import sys, os, re, shutil, Bio, csv
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
				GC3 = splits[5]
				Aratio = splits[10]

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


def CAI_a_content(HE_genes, CAIvals, geneSeq):


	fourthA = {}
	fourthNotA = {}

	meanFourthA = {}
	meanFourthNotA = {}

	# For each of the genomes with CAI values
	for acc in CAIvals:

		fourthA[acc] = []
		fourthNotA[acc] = []

		# For each of the gene loci in the genome
		for locus in CAIvals[acc]:

			# Ensure the gnee is not one of the reference genes
			if locus not in HE_genes[acc]:


				if CAIvals[acc][locus] != 0:

					# If that locus is also in the gene sequences
					if locus in geneSeq[acc]:

						CAIval = CAIvals[acc][locus]
						seq = geneSeq[acc][locus]

						if seq[3] == "A":
							fourthA[acc].append(CAIval)
						else:
							fourthNotA[acc].append(CAIval)

		if fourthA[acc]:
			meanA = sum(fourthA[acc]) / float(len(fourthA[acc]))
			meanFourthA[acc] = meanA

		if fourthNotA[acc]:
			meanNotA = sum(fourthNotA[acc]) / float(len(fourthNotA[acc]))
			meanFourthNotA[acc] = meanNotA



	return meanFourthA, meanFourthNotA



def write_CAI_a_not_a(GC3s, meanA, meanNotA, CAIaNotaPath, CAIaNotaPathRestrictGC3):

	file = open(CAIaNotaPath, "w")
	file.write("acc,gc3,meanCAI_a,meanCAI_not_a\n")

	# file2 = open(CAIaNotaPathRestrictGC3, "w")
	# file2.write("acc,gc3,meanCAI_a_restricted,meanCAI_not_a_restricted\n")


	for acc in meanA:

		if acc in GC3s:

			fileLine = "%s,%s,%s,%s\n" % (acc,GC3s[acc], meanA[acc], meanNotA[acc])
			file.write(fileLine)

			# # Restrict the genomes without extreme GC3
			# if float(GC3s[acc]) > 0.2 and float(GC3s[acc]) < 0.9:
			#
			# 	fileLine2 = "%s,%s,%s,%s\n" % (acc,GC3s[acc], meanA[acc], meanNotA[acc])
			# 	file2.write(fileLine2)



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


# Compare the A content between high and low CAI genes
def compare_high_low_CAI(HE_genes, CAIvals, geneSeq, GC3s, highLowPath, highLowPathRestrictedGC3):




	# Get the max/main values of the highly expressed genes
	maxMinCAIs = {}
	for acc in HE_genes:




		maxMinCAIs[acc] = [0,1]


		for locus in HE_genes[acc]:



			if locus in CAIvals[acc]:

				# print acc, locus

				CAIval = CAIvals[acc][locus]




				if CAIval > maxMinCAIs[acc][0]:
					maxMinCAIs[acc][0] = CAIval
				if CAIval < maxMinCAIs[acc][1]:
					maxMinCAIs[acc][1] = CAIval

		# print maxMinCAIs[acc][0], maxMinCAIs[acc][1]

	# print len(HEgenes)
	# for acc in maxMinCAIs:
	# 	print acc, maxMinCAIs[acc][0], maxMinCAIs[acc][1]
	# print len(maxMinCAIs)

	highCAIs = {}
	otherCAIs = {}



	# Get the genes with a CAI in the range of HE genes
	for acc in CAIvals:

		highCAIs[acc] = []
		otherCAIs[acc] = []



		for locus in CAIvals[acc]:


			# Ensure the locus is not one of the reference genes
			if locus not in HE_genes[acc]:


				if locus in CAIvals[acc]:



					if CAIvals[acc][locus] <= maxMinCAIs[acc][0] and CAIvals[acc][locus] >= maxMinCAIs[acc][1]:
						highCAIs[acc].append(locus)
					else:
						otherCAIs[acc].append(locus)

		# print highCAIs[acc]
		# print otherCAIs[acc]

	# for acc in highCAIs:
	# 	print highCAIs[acc], otherCAIs[acc]
	# print len(highCAIs)


	# Get the fourth site of each of the high CAI genes and other
	highCAIfourths = {}
	otherCAIfourths = {}

	for acc in highCAIs:

		highCAIfourths[acc] = []
		for locus in highCAIs[acc]:

			if locus in geneSeq[acc]:
				seq = geneSeq[acc][locus]
				highCAIfourths[acc].append(seq[3])

	for acc in otherCAIs:

		# print acc
		otherCAIfourths[acc] = []


		for locus in otherCAIs[acc]:



			# print locus

			if locus in geneSeq[acc]:
				seq = geneSeq[acc][locus]
				otherCAIfourths[acc].append(seq[3])


	highCAIprops = {}
	otherCAIprops = {}

	for acc in highCAIfourths:



		if len(highCAIfourths[acc]) != 0:



			highPropA = highCAIfourths[acc].count("A") / float(len(highCAIfourths[acc]))
		else:
			highPropA = 0

		highCAIprops[acc] = highPropA


	for acc in otherCAIfourths:
		if len(otherCAIfourths[acc]) != 0:
			otherPropA = otherCAIfourths[acc].count("A") / float(len(otherCAIfourths[acc]))
		else:
			otherPropA = 0


		otherCAIprops[acc] = otherPropA


	highOtherFile = open(highLowPath, "w")
	highOtherFile.write("acc,gc3,prop_a_high,prop_a_other\n")

	# highOtherFile2 = open(highLowPathRestrictedGC3, "w")
	# highOtherFile2.write("acc,prop_a_high_restricted,prop_a_other_restricted\n")

	for acc in highCAIprops:


		if acc in GC3s:

			# print acc, highCAIprops[acc], otherCAIprops[acc]

			if highCAIprops[acc] != 0 and  otherCAIprops[acc] != 0:



				fileLine = "%s,%s,%s,%s\n" % (acc, GC3s[acc], highCAIprops[acc], otherCAIprops[acc])
				highOtherFile.write(fileLine)

				# if float(GC3s[acc]) > 0.2 and float(GC3s[acc]) < 0.9:
				# 	fileLine2= "%s,%s,%s\n" % (acc, highCAIprops[acc], otherCAIprops[acc])
				# 	highOtherFile2.write(fileLine2)



	highOtherFile.close()


	return maxMinCAIs








def compare_codonw(CAIvals, comparePath):

	codonWCAIs = {}

	codonWpath = expression_genomes_dir + "AE005174/AE005174codonw.txt"
	with open(codonWpath, "U") as myfile:

		for line in myfile:

			atts = line.split("\t")
			gene = atts[0]
			CAI = atts[8]
			codonWCAIs[gene] = CAI


	compareFile = open(comparePath, "w")
	compareFile.write("gene,my_cai,codonw_cai\n")
	for gene in CAIvals["AE005174"]:

		fileLine = "%s,%s,%s\n" % (gene, CAIvals["AE005174"][gene], codonWCAIs[gene])
		compareFile.write(fileLine)



def compare_a_content_low_mean_cai(CAIvals, geneSeq, GC3s, HE_genes, meanCAIs, low_mean_cai_a_content_path, high_mean_cai_a_content_path):


	highCAIsHigh = {}
	otherCAIsHigh = {}
	highCAIsLow = {}
	otherCAIsLow = {}

	for acc in meanCAIs:

		# If the acc mean CAI is less than 0.6, continue
		if meanCAIs[acc] > 0.6:

			highCAIsHigh[acc] = []
			otherCAIsHigh[acc] = []

			# Get the CAI values for the highly expressed genes in the genomes with a
			# mean CAI < 0.6
			for gene in HE_genes[acc]:
				if gene in CAIvals[acc]:
					highCAIsHigh[acc].append(float(CAIvals[acc][gene]))

		else:

			highCAIsLow[acc] = []
			otherCAIsLow[acc] = []

			for gene in HE_genes[acc]:
				if gene in CAIvals[acc]:
					highCAIsLow[acc].append(float(CAIvals[acc][gene]))


	# Get the minimum CAI value for the highly expressed genes
	# in both high CAI mean genomes and others
	minHighCAI = {}
	for acc in highCAIsHigh:
		if highCAIsHigh[acc]:
			minHighCAI[acc] = min(highCAIsHigh[acc])

	minLowCAI = {}
	for acc in highCAIsLow:
		if highCAIsLow[acc]:
			minLowCAI[acc] = min(highCAIsLow[acc])



	fourthHighCAIhigh = {}
	fourthOtherCAIhigh = {}

	# Get the fourth positions of genes in the high mean CAI group
	for acc in minHighCAI:

		fourthHighCAIhigh[acc] = []
		fourthOtherCAIhigh[acc] = []

		for gene in CAIvals[acc]:


			# Ensure the gene is not one of the reference genes
			if gene not in HE_genes[acc]:


				if float(CAIvals[acc][gene]) >= minHighCAI[acc]:

				# 	print "true", CAIvals[acc][gene], minHighCAI[acc]

				# else:
				# 	print "false", CAIvals[acc][gene], minHighCAI[acc]

					if gene in geneSeq[acc]:
						fourthHighCAIhigh[acc].append(geneSeq[acc][gene][3])
				else:
					if gene in geneSeq[acc]:
						fourthOtherCAIhigh[acc].append(geneSeq[acc][gene][3])



	fourthHighCAIlow = {}
	fourthOtherCAIlow = {}
	# Get the fourth positions of genes in the low mean CAI group

	for acc in minLowCAI:

		fourthHighCAIlow[acc] = []
		fourthOtherCAIlow[acc] = []

		for gene in CAIvals[acc]:


			# Ensure the gene is not one of the reference genes
			if gene not in HE_genes[acc]:
				if float(CAIvals[acc][gene]) >= minLowCAI[acc]:

					if gene in geneSeq[acc]:
						fourthHighCAIlow[acc].append(geneSeq[acc][gene][3])
				else:
					if gene in geneSeq[acc]:
						fourthOtherCAIlow[acc].append(geneSeq[acc][gene][3])




	file = open(high_mean_cai_a_content_path, "w")
	file.write("acc,gc3,mean_cai,prop_a_high,prop_a_other\n")
	for acc in fourthHighCAIhigh:

		if acc in GC3s:

			fileLine = "%s,%s,%s,%s,%s\n" % (acc, GC3s[acc], meanCAIs[acc], fourthHighCAIhigh[acc].count("A") / float(len(fourthHighCAIhigh[acc])), fourthOtherCAIhigh[acc].count("A") / float(len(fourthOtherCAIhigh[acc])))
			file.write(fileLine)
	file.close()



	file2 = open(low_mean_cai_a_content_path, "w")
	file2.write("acc,gc3,mean_cai,prop_a_high,prop_a_other\n")
	for acc in fourthHighCAIlow:

		if acc in GC3s:

			fileLine = "%s,%s,%s,%s,%s\n" % (acc, GC3s[acc], meanCAIs[acc], fourthHighCAIlow[acc].count("A") / float(len(fourthHighCAIlow[acc])), fourthOtherCAIlow[acc].count("A") / float(len(fourthOtherCAIlow[acc])))
			file2.write(fileLine)
	file2.close()

#################################

def main():


	fileNumber = 0

	# Get a list of the genomes containing CAI information
	genomes, numGenomes = getGenomes(expression_genomes_dir)



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

			print "-" * 10
			print acc


			# Get the CAI values for each genome
			getCAIs(acc, CAIvals)
			geneSeq = get_geneSeq(acc, geneSeq)





	# Compare the CAI values of those calculated using HE genes
	# and the default values given by CodonW
	print "\n\nGetting the CAI values for E. coli\n\n"
	comparePath = expression_dir + "_ecoli_codonw_compare.csv"
	compare_codonw(CAIvals, comparePath)


	# Calculate the mean CAI for each genome
	print "\n\nCalculating the mean CAI for each genome\n\n"

	fileNumber = 0

	# Get a list of the highly expressed ribosomal genes used as
	# a reference set
	HE_genes = get_HE_genes(CAIvals)

	# Calculate the mean CAI for the genome
	meanCAIs = calcMeanCAIs(HE_genes, CAIvals, acc)


	CAIGCpath = expression_dir + "_CAI_a_ratios.csv"
	write_mean_CAIGC(CAIGCpath, meanCAIs, GC3s, Aratios)

	# Calculate the mean CAI for each genome
	print "\n\nCalculating CAIs for genes with and without +4 A\n\n"
	meanA, meanNotA = CAI_a_content(HE_genes, CAIvals, geneSeq)
	CAIaNotaPath = expression_dir + "_CAI_a_notA.csv"
	CAIaNotaPathRestrictGC3 = expression_dir + "_CAI_a_nota_restrictedGC3.csv"
	write_CAI_a_not_a(GC3s, meanA, meanNotA, CAIaNotaPath, CAIaNotaPathRestrictGC3)


	# Compare the highly expressed genes with others
	print "\n\nComparing A content between high CAI genes and low CAI genes\n\n"
	highLowPath = expression_dir + "_high_low_CAI_compare.csv"
	highLowPathRestrictedGC3 = expression_dir + "high_low_CAI_compare_restrictedGC3.csv"
	maxMinCAIs = compare_high_low_CAI(HE_genes, CAIvals, geneSeq, GC3s, highLowPath, highLowPathRestrictedGC3)


	# Compare the proportion of +4A in genes for genomes with high/low mean CAI
	print "\n\nComparing A content between high CAI genes and low CAI genes for genomes with high/low mean CAI\n\n"
	low_mean_cai_a_content_path = expression_dir + "_lowMeanCAIcompareFourthA.csv"
	high_mean_cai_a_content_path = expression_dir + "_highMeanCAIcompareFourthA.csv"
	compare_a_content_low_mean_cai(CAIvals, geneSeq, GC3s, HE_genes, meanCAIs, low_mean_cai_a_content_path, high_mean_cai_a_content_path)




#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
