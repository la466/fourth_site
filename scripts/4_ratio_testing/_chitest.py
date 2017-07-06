#!/usr/bin/python

# Script number: 			4.2
# File: 					2 of 6
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 4.1
# Prerequisite file(s):		_site_SITE_ratios.csv, site_SITE_ratios_t4.csv
# Input(s):					-site SITE_TO_CALULATE(int)
# Output(s):				_chisquare_site_SITE.csv, _chisquare_site_SITE_t4.csv

import sys, os, csv
import rpy2.robjects as robjects
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description = "Help for Chi square calculation script.", formatter_class=RawTextHelpFormatter)
parser.add_argument("site", help="Enter the site you wish to calculate chi square values for.", type=int)
args = parser.parse_args()


if args.site:

	if float(args.site) > 12:
		print "\nWARNING: You have chosen a site greater than 12. Please re-run script 4.2 (_chitest.py) using another site.\n"
		sys.exit()
	else:
		print "\nCalculating chi square values for site %d\n" % args.site
else:
	sys.exit()



##########################
# VARIABLES
##########################



summaryFile = "outputs/ratio_testing/site_%d/_site_%s_ratios.csv" % (args.site, args.site)
ChiSquareFilePath = "outputs/ratio_testing/site_%d/_chisquare_site_%s.csv" % (args.site, args.site)

summaryFileT4 = "outputs/ratio_testing/site_%d/_site_%s_ratios_t4.csv" % (args.site, args.site)
ChiSquareFilePathT4 = "outputs/ratio_testing/site_%d/_chisquare_site_%s_t4.csv" % (args.site, args.site)


if not os.path.exists(summaryFile):
	print "\nWARNING: Cannot find the ratio file '%s' used to perform chi square test. Please re-run script 4.1 (_calculate_ratios.py) using -cds_site %d.\n" % (summaryFile, args.site)
	sys.exit()

##########################
# FUNCTIONS
##########################

# Calculate chiqaure value
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
			gc = row[4]
			gc3 = row[5]

			a_ratio = row[10]

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
	chisqFileT4 = open(ChiSquareFilePathT4, "w")

	fileHead = "genus,acc_num,num_genes,gc,gc3,a_ratio,obs_a4,exp_a4,oba_not_a4,exp_not_a4,chival,pval\n"

	chisqFile.write(fileHead)
	chisqFileT4.write(fileHead)

	chiSquareTest(summaryFile,chisqFile)
	chiSquareTest(summaryFileT4,chisqFileT4)

	chisqFile.close()
	chisqFileT4.close()


#################################

if __name__ == "__main__":
	main()
