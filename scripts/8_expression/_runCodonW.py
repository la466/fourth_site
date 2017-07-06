#!/usr/bin/python

# Script number: 			8.3
# File: 					3 of 5
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 8.1, 8.2
# Prerequisite file(s):		ACC/ACC.txt, ACC/ACC_he.txt
# Description: 				Run CodonW on all genomes. Use highly expressed genes to get indicies for the rest of the genes. For E.coli, also use CodonW indices.
# Output file(s):			ACC/ACC.txt, ACC/ACC_he.txt



import sys, os, re, shutil, Bio, csv
from sys import argv


##########################
# VARIABLES
##########################


expression_genomes_dir = "outputs/expression/genomes/"



testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE005174"]

coaFiles = ["cai.coa", "cbi.coa", "coa_raw", "codon.coa", "cusort.coa", "eigen.coa", "fop.coa", "genes.coa", "hilo.coa", "summary.coa"]



##########################
# FUNCTIONS
##########################


# Create a list of all the good files
def getGenomes(direc):

	files = []

	# For each file in the downloaded genomes
	for eachfile in os.listdir(direc):

		if eachfile != '.DS_Store':
			# if eachfile in testFile:
				files.append(eachfile)

	return files



def CAIanalysis(fileNumber, genome,numGenomes):

	print '-' * 10
	print "%s of %s" % (fileNumber, numGenomes)
	print genome

	HEfile = expression_genomes_dir + genome + "/" + genome + "_he.txt"
	allGenesFile = expression_genomes_dir + genome + "/" + genome + ".txt"

	codonWcommand = "packages/codonW/./codonw %s -coa_cu -coa_num %s -nomenu -silent" % (HEfile, "100%")
	os.system(codonWcommand)

	coaPath = expression_genomes_dir + genome + "/" + "coa/"

	if not os.path.exists(coaPath):
		os.makedirs(coaPath)

	for file in coaFiles:
		filePath = "bacterial_genomes/" + file

		coaFile = coaPath + file
		if os.path.exists(coaFile):
			os.remove(coaFile)


		shutil.move(file, coaPath)

	codonWcommand2 = "packages/codonW/./codonw %s -all_indices -fop_file %s -cai_file %s -cbi_file %s -nomenu -silent" % (allGenesFile, coaPath + "fop.coa", coaPath + "cai.coa", coaPath + "cbi.coa")
	# codonWcommand2 = "codonW/./codonw %s -all_indices -cai_file %s -nomenu -silent" % (allGenesFile)

	os.system(codonWcommand2)


	if genome == "AE005174":
		codonWcommand3 = "packages/codonW/./codonw %s %s -all_indices -nomenu -silent" % (allGenesFile, expression_genomes_dir + genome + "/" + genome + "codonw.txt")
		os.system(codonWcommand3)

#################################

def main():


	fileNumber = 0

	genomes = getGenomes(expression_genomes_dir)


	for genome in genomes:
		fileNumber +=1
		CAIanalysis(fileNumber, genome, len(genomes))


#################################


if __name__ == "__main__":
	main()
