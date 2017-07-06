#!/usr/bin/python

# Script number: 			4.4
# File: 					4 of 6
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 4.3
# Prerequisite file(s):		_site_4_no_overlap.csv
# Input(s):
# Output(s):				_chisquare_site_4_no_overlap.csv

import sys, os, csv
import rpy2.robjects as robjects


print "\nPerforming Chi Square test on the genome ratios with no 4 site overlaps.\n"


##########################
# VARIABLES
##########################


summaryFile = "outputs/ratio_testing/site_4/_site_4_no_overlap.csv"
ChiSquareFilePath = "outputs/ratio_testing/site_4/_chisquare_site_4_no_overlap.csv"


##########################
# FUNCTIONS
##########################

def chisquare(obs,exp):
	chisq = ((obs-exp)**2)/float(exp)
	return chisq


def chiSquareTest(filePath,outputFile):

	rawFile = csv.reader(open(filePath, "rU"), dialect=csv.excel_tab)

	rowNumber = 0

	for row in rawFile:

		row = row[0].split(',')

		if rowNumber > 0:
			genus = row[0]
			acc_num = row[1]
			num_genes = int(row[3])
			gc3 = row[4]

			a_ratio = row[9]


			obs_a4 = int(row[6])
			obs_prop_a = float(row[8])
			exp_a4 = obs_prop_a*num_genes

			obs_c4 = int(row[21])
			obs_prop_c = float(row[23])
			exp_c4 = obs_prop_c*num_genes

			obs_t4 = int(row[11])
			obs_prop_t = float(row[13])
			exp_t4 = obs_prop_t*num_genes

			obs_g4 = int(row[16])
			obs_prop_g = float(row[18])
			exp_g4 = obs_prop_g*num_genes

			obs_not_a = obs_c4 + obs_g4 + obs_t4
			exp_not_a = exp_c4 + exp_g4 + exp_t4


			chival = (((obs_a4-exp_a4)**2)/exp_a4) + (((obs_not_a-exp_not_a)**2)/exp_not_a)
			pval = robjects.r('pchisq(%s,1,lower.tail=FALSE)' % chival)[0]



			# pvalA = robjects.r('pchisq(%s,1,lower.tail=FALSE)' % chisqA)[0]
			# pvalT = robjects.r('pchisq(%s,1,lower.tail=FALSE)' % chisqT)[0]
			# pvalG = robjects.r('pchisq(%s,1,lower.tail=FALSE)' % chisqG)[0]
			# pvalC = robjects.r('pchisq(%s,1,lower.tail=FALSE)' % chisqC)[0]
			# totalpval = robjects.r('pchisq(%s,3,lower.tail=FALSE)' % totalChisq)[0]

			outputLine = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (genus,acc_num,num_genes,gc3,a_ratio,obs_a4,exp_a4,obs_not_a,exp_not_a,chival,pval)
			outputFile.write(outputLine)


		rowNumber += 1


def main():
	chisqFile = open(ChiSquareFilePath, "w")

	fileHead = "genus,acc_num,num_genes,gc3,a_ratio,obs_a4,exp_a4,oba_not_a4,exp_not_a4,chival,pval\n"

	chisqFile.write(fileHead)

	chiSquareTest(summaryFile,chisqFile)

	chisqFile.close()


#################################


if __name__ == "__main__":
	main()
