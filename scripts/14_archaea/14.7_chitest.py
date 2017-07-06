#!/usr/bin/python

# Script number: 			14.7
# File: 					7 of 7
# Prerequisite script(s):	14.1, 14.2, 14.3, 14.4, 14.5, 14.6
# Prerequisite file(s):		_site_SITE_ratios.csv, site_SITE_ratios_t4.csv
# Input(s):					-site SITE_TO_CALULATE(int)
# Output(s):				_chisquare_site_SITE.csv, _chisquare_site_SITE_t4.csv


import sys, os, csv
import rpy2.robjects as robjects
import argparse
from argparse import RawTextHelpFormatter




##########################
# VARIABLES
##########################



summaryFile = "outputs/archaea/ratio_testing/site_4/_site_4_ratios_archaea.csv"
ChiSquareFilePath = "outputs/archaea/ratio_testing/site_4/_chisquare_site_4_archaea.csv"





##########################
# FUNCTIONS
##########################

# Calculate chiqaure value
def chisquare(obs,exp):

	if exp != 0:
		chisq = ((obs-exp)**2)/float(exp)
	else:
		chisq = 0

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
			gc = row[4]
			gc3 = row[5]

			a_ratio = row[10]
			# c_ratio = row[26]
			# t_ratio = row[16]
			# g_ratio = row[21]

			obs_a4 = int(row[7])
			obs_prop_a = float(row[9])
			exp_a4 = obs_prop_a*num_genes

			obs_c4 = int(row[22])
			obs_prop_c = float(row[24])
			exp_c4 = obs_prop_c*num_genes

			obs_t4 = int(row[12])
			obs_prop_t = float(row[14])
			exp_t4 = obs_prop_t*num_genes

			obs_g4 = int(row[17])
			obs_prop_g = float(row[19])
			exp_g4 = obs_prop_g*num_genes

			obs_not_a = obs_c4 + obs_g4 + obs_t4
			exp_not_a = exp_c4 + exp_g4 + exp_t4


			chival = (((obs_a4-exp_a4)**2)/exp_a4) + (((obs_not_a-exp_not_a)**2)/exp_not_a)
			pval = robjects.r('pchisq(%s,1,lower.tail=FALSE)' % chival)[0]


			outputLine = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (genus,acc_num,num_genes,gc,gc3,a_ratio,obs_a4,exp_a4,obs_not_a,exp_not_a,chival,pval)
			outputFile.write(outputLine)


		rowNumber += 1


def main():
	chisqFile = open(ChiSquareFilePath, "w")


	fileHead = "genus,acc_num,num_genes,gc,gc3,a_ratio,obs_a4,exp_a4,oba_not_a4,exp_not_a4,chival,pval\n"
	chisqFile.write(fileHead)


	chiSquareTest(summaryFile,chisqFile)


	chisqFile.close()



#################################

if __name__ == "__main__":
	main()
