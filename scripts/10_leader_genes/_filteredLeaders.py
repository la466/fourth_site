#!/usr/bin/python

# Script number: 			10.2
# File: 					4 of 4
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 10.1
# Prerequisite file(s):		_leaderGenes.csv, _filteredGenes.csv, table4genomes.txt, _site_4_ratios.csv
# Description: 				Extract leader genes based upon previous CDS filtering restrictions
# Output file(s):			_leaderGenes.csv, _number_leaders.csv



import re

from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate



##########################
# VARIABLES
##########################


leadersFile = "outputs/leader_genes/_leaderGenes.csv"
filteredCDSfiles = "outputs/gene_filtering/_filteredGenes.csv"

filteredLeadersFile = "outputs/leader_genes/_leader_genes_filtered.csv"
numberLeadersFile = "outputs/leader_genes/_number_leaders.csv"
fourthSites = "outputs/ratio_testing/site_4/_site_4_ratios.csv"

table4_genomes_path = "outputs/gene_filtering/table4genomes.txt"



##########################
# FUNCTIONS
##########################


def get_table4_genomes(file_path):

	with open(file_path, "U") as myfile:

		table4_genomes = []

		lines = myfile.readlines()
		for line in lines:
			line = line.replace("\n", "")
			table4_acc = line.split(",")[1]

			table4_genomes.append(table4_acc)

		return table4_genomes


def get_GC(fourthSites):

	GC3s = {}
	GCs = {}
	rowNum = 0

	with open(fourthSites) as myfile:

		for row in myfile:

			rowNum += 1

			if rowNum > 1:

				splits = row.split(",")

				acc = splits[1]
				gc3 = splits[5]
				gc = splits[4]

				GC3s[acc] = gc3
				GCs[acc] = gc


	return GC3s, GCs


# Get the list of filtered genes
def get_filtered_cds(file, table4_genomes):

	filteredGenes = {}
	rowNum = 0

	# Open the file containing the filtered genes
	with open(file, "U") as myfile:

		# For each row in the file
		for row in myfile:

			rowNum += 1

			# Exclude the header row
			if rowNum > 1:

				# Split the row
				genome = row.split(',')

				# Get the accession number of the genome
				acc = genome[1]

				if acc not in table4_genomes:

					filteredGenes[acc] = []

					# For each of the filtered gene in the genome
					for gene in range(2, len(genome)):

						filtered_locus = genome[gene].strip('\n')

						# Append the filtered genes to the dictionary
						filteredGenes[acc].append(filtered_locus)


	return filteredGenes


# Determine if the CDS for the leader gene is in the filtered genes
def filter_leader_genes(filteredCDS, leadersFile, filteredLeadersFile, numberLeadersFile, gc3s, gcs, table4_genomes):

	totalLeaders = {}

	filterFile = open(filteredLeadersFile, 'w')
	filterFile.write("acc,locus,cds,leader\n")

	numberFile = open(numberLeadersFile, 'w')
	#numberFile.write("acc,phylum,leaders,total,prop,gc3,gc\n")
	numberFile.write("acc,leaders,total,prop,gc3,gc\n")

	rowNum = 0

	# For each of the identified leader genes
	with open(leadersFile, "U") as myfile:

		for row in myfile:

			rowNum +=1

			# Ignore the header line
			if rowNum > 1:

				leader = row.split(',')
				acc = leader[0]

				if acc not in table4_genomes and acc in filteredCDS:

					locus = leader[1]
					cds = leader[2]
					leaderregion = leader[3].strip("\n")

					# If the CDS locus is one of the filtered genes, keep the leader gene
					if locus in filteredCDS[acc]:

						fileLine  = "%s,%s,%s,%s\n" % (acc,locus,cds,leaderregion)
						filterFile.write(fileLine)

						if acc in totalLeaders:
							totalLeaders[acc] += 1
						else:
							totalLeaders[acc] = 1


	print "%s genomes\n" % len(totalLeaders)
	genome_count = 0

	for acc in totalLeaders:


		genome_count +=1
		print acc
		print "%s of %s" % (genome_count, len(totalLeaders))

		leaders = totalLeaders[acc]
		cds = len(filteredCDS[acc])
		prop = leaders / float(cds)

		#fileLine = "%s,%s,%s,%s,%s,%s,%s\n" % (acc, phylum, leaders, cds, prop, gc3s[acc], gcs[acc])
		fileLine = "%s,%s,%s,%s,%s,%s\n" % (acc, leaders, cds, prop, gc3s[acc], gcs[acc])

		numberFile.write(fileLine)



	filterFile.close()
	numberFile.close()

#################################

def main():

	table4_genomes = get_table4_genomes(table4_genomes_path)

	gc3s, gcs = get_GC(fourthSites)

	print "\nFiltering the CDSs with a potential leader gene\n"
	filteredGenes = get_filtered_cds(filteredCDSfiles, table4_genomes)

	filter_leader_genes(filteredGenes, leadersFile, filteredLeadersFile, numberLeadersFile, gc3s, gcs, table4_genomes)




#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
