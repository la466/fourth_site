#!/usr/bin/python

# Script number: 			5.3
# File: 					3 of 3
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 5.1, 5.2
# Prerequisite file(s):		_codon_gc_proportions.csv
# Description: 				Calculate the variance in GC content for each position in codon 2-30.
# Output file(s):			_codon_gc_variance.csv

import sys, os, csv, re, shutil
from sys import argv
from collections import deque


print "\nCalculating the GC variance for each position in codons 2-30.\n"

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
	nucleotidePath = baseCompDir + "_codon_gc_variance.csv"
	nucleotideFile = open(nucleotidePath, "w")

	return nucleotidePath, nucleotideFile



def fileHead(file):

	codonLine = ""
	codonLine += "codon,Base 1,Base 2,Base 3\n"
	file.write(codonLine)




def proportion(number, total):

	if total != 0:
		pc = number / float(total)
	else:
		pc = 0
	return pc





# Get all the good genes
def baseUsage(nucleotideFile):



	

	nucleotidesCompFile = baseCompDir + "/_codon_gc_proportions.csv"
	with open(nucleotidesCompFile, 'U') as myfile:
		# read = csv.reader(myfile)
		

		codons = {}
		for codon in codonRange:
			codons[codon] = {}
			for base in basePos:
				codons[codon][base] = 0


		totals = deque(csv.reader(myfile), 1)[0]
		totals = totals[4:]
		
		i = 0

		codon = 0

		for i in range(0,len(totals)):
			
			if (i+3) % 3 == 0:
				z=1
				codon = (i+6)/3
				codons[codon][z] = totals[i]
				z += 1
			elif (i+2) % 3 == 0:
				codons[codon][z] = totals[i]
				z+=1
			elif (i+1) % 3 == 0:				
				codons[codon][z] = totals[i]

			
	for codon in codons:
		codonLine = ""	
		codonLine += "%s,%s,%s,%s\n"  % (codon,codons[codon][1], codons[codon][2], codons[codon][3])
		nucleotideFile.write(codonLine)	

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

