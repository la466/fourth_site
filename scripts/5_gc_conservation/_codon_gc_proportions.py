#!/usr/bin/python

# Script number: 			5.2
# File: 					2 of 3
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 5.1
# Prerequisite file(s):		_codon_nucleotide_proportions.csv
# Description: 				Calculate the GC content for each position in codon 2-30.
#							Only use genes that are 30+ codons long
# Output file(s):			_codon_gc_proportions.csv

import sys, os, csv, re, shutil
from sys import argv
import rpy2.robjects as robjects


##########################
# VARIABLES
##########################

# Set the directory containing the sorted and parsed genomes
baseCompDir = 'outputs/nucleotide_conservation/'

testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE000511"]

# Define the bases
bases = ["A", "C", "T", "G"]

codonRange = range(2,31)
basePos = range(1,4)


##########################
# FUNCTIONS
##########################

# Define and open the output files
def newFile():

	# Output file with genes passing the filter
	nucleotidePath = baseCompDir + "_codon_gc_proportions.csv"
	nucleotideFile = open(nucleotidePath, "w")

	return nucleotidePath, nucleotideFile



def fileHead(file):

	codonLine = ""
	codonLine += "genus,acc,num_genes,genome_gc"

	for codon in codonRange:
		for pos in basePos:
			codonLine += ",codon%s_base%s" % (codon, pos)
	
	codonLine += "\n"
	file.write(codonLine)




def proportion(number, total):

	if total != 0:
		pc = number / float(total)
	else:
		pc = 0
	return pc





# Get all the good genes
def baseUsage(nucleotideFile):



	goodGenes = {}

	nucleotidesCompFile = baseCompDir + "_codon_nucleotide_proportions.csv"
	with open(nucleotidesCompFile, 'U') as myfile:
		read = csv.reader(myfile)
		
		GCusage = {}
		rowNum = 0

		totalUsage = {}


		for codon in codonRange:
			GCusage[codon] = {}
			totalUsage[codon] = {}
			for base in basePos:
				GCusage[codon][base] = 0
				totalUsage[codon][base] = []

		

		for row in read:		

			if rowNum > 0:

				accLine = ""
				
				genus = row[0]
				acc = row[1]
				numGenes = row[2]
				genomeGC = row[3]

				print "-" * 20
				print acc
				print "Genome %s" % (rowNum)

				
				accLine += "%s,%s,%s,%s" % (genus,acc,numGenes,genomeGC)

				
				for i in codonRange:
					base1C = row[((i-1)*12)-7]
					base1G = row[((i-1)*12)-5]
					base2C = row[((i-1)*12)-3]
					base2G = row[((i-1)*12)-1]
					base3C = row[((i-1)*12)+1]
					base3G = row[((i-1)*12+3)]


					GCusage[i][1] = float(base1G) + float(base1C)
					GCusage[i][2] = float(base2G) + float(base2C)
					GCusage[i][3] = float(base3G) + float(base3C)

					# Add the GC content for that codons position to the list
					totalUsage[i][1].append(float(base1G) + float(base1C))
					totalUsage[i][2].append(float(base2G) + float(base2C))
					totalUsage[i][3].append(float(base3G) + float(base3C))



				for codon in GCusage:
					for position in GCusage[codon]:
						accLine += ",%s" % GCusage[codon][position]

				accLine += "\n"

				nucleotideFile.write(accLine)



						# print v
			rowNum += 1
			
		GCvars = {}
		for codon in codonRange:
			GCvars[codon] = {}
			for position in basePos:
				GCvars[codon][position] = 0


		for codon in totalUsage:
			for position in totalUsage[codon]:

				# print codon, position, totalUsage[codon][position]
				test = totalUsage[codon][position]
				
				# Get the codon position GC usages as an R vector
				v = robjects.FloatVector(test)
				

				rvar = robjects.r['var']
				
				# Calculate the variance for each position in each codon
				var = rvar(v)[0]

				GCvars[codon][position] = var
		

		totalLine = "Total,-,-,-"

		for codon in GCvars:
			for position in GCvars[codon]:
				# print codon, position, GCvars[codon][position]
				totalLine += ",%s" % GCvars[codon][position]

		nucleotideFile.write(totalLine)
#################################

def main():

	# Set the number of files to 0
	fileNumber = 0
	

	
	# Open the new files
	nucleotidePath, nucleotideFile = newFile()

	# Write the heading for the new files
	fileHead(nucleotideFile)

	baseUsage(nucleotideFile)


	# Close both of the output files
	nucleotideFile.close()

#################################

if __name__ == "__main__":
	main()

