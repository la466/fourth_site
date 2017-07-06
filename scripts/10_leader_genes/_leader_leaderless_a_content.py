#!/usr/bin/python

# Script number: 			10.4
# File: 					4 of 4
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 10.1, 10.2
# Prerequisite file(s):		_leader_genes_filtered.csv, _filteredGenes.csv
# Description: 				Look at the A content of CDSs with and without leader genes
# Output file(s):			_leader_leaderless_a_content.csv


import re

from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate

import re



##########################
# VARIABLES
##########################

filteredLeadersFile = "outputs/leader_genes/_leader_genes_filtered.csv"
emblPath = "bacterial_genomes/good_genomes/"
leaderFilePath = "outputs/leader_genes/_leader_leaderless_a_content.csv"
table4_genomes_path = "outputs/gene_filtering/table4genomes.txt"


# testFile = ["AE005174", "AE000511", "AE000512"]


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


# Get the leader genes from the file
def get_leaders(file, table4_genomes):

	leaders = {}
	rowNum = 0

	# Open the file containing the filtered genes
	with open(file, "U") as myfile:

		# For each row in the file
		for row in myfile:

			rowNum += 1

			# Exclude the header row
			if rowNum > 1:

				# Split the row
				leadergene = row.split(',')

				# Get the accession number of the genome
				acc = leadergene[0]

				if acc not in table4_genomes:
					locus = leadergene[1]
					cds = leadergene[2]
					leaderRegion = leadergene[3].strip("\n")

					if acc in leaders:
						leaders[acc][locus] = [cds,leaderRegion]
					else:
						leaders[acc] = {}
						leaders[acc][locus] = [cds,leaderRegion]

	return leaders



def getFilteredGenes():

	goodGenes = {}

	filteredGenesFile = "outputs/gene_filtering/_filteredGenes.csv"
	with open(filteredGenesFile, 'U') as myfile:


		genomes = {}

		rowNum = 0

		# For each row in the file
		for row in myfile:
			rowNum += 1

			# Ignore the header
			if rowNum > 1:


				# Split the row
				splits = row.split(",")

				# Get the accession number
				acc = splits[1]

				# Get each of the filtered genes
				genes = []
				for i in range(2, len(splits)):

					locus = splits[i]
					if i == len(splits)-1:
						locus = locus.strip("\n")

					genes.append(locus)

					# print genes

				genomes[acc] = genes

		return genomes



def get_filtered_cds(leaders, filteredCDS):


	genomeNo = 0
	numGenomes = len(leaders)
	cds = {}

	for acc in leaders:

		# if acc in testFile:
			genomeNo += 1
			print acc
			print "%s of %s" % (genomeNo, numGenomes)
			print "-" * 20


			# Get the path to the embl file
			emblFile = emblPath + acc + ".embl"

			# Get the genomes raw sequence
			for seq_record in SeqIO.parse(emblFile, "embl"):
				raw_sequence = seq_record.seq


			# Read the file using Biopython
			record = SeqIO.read(emblFile, "embl")



			locusTags = []

			cdsNum = 0

			cds[acc] = {}

			# Retrieve the CDS for each genomes
			for feature in range(1,len(record.features)):

				# If the feature type is CDS
				if record.features[feature].type == "CDS":

					cdsNum += 1

					# Get the locus tag
					locus_tag = record.features[feature].qualifiers.get('locus_tag')
					translation_table = record.features[feature].qualifiers.get('transl_table')



					if not locus_tag:
						locus_tag = "no_tag_%s" % cdsNum
					else:
						locus_tag = locus_tag[0]


					samelocus = 0

					for singleLocus in locusTags:
						if locus_tag in singleLocus:
							samelocus += 1

					if samelocus > 0:
						locus_tag = "%s_%s" % (locus_tag, samelocus + 1)


					# print locus_tag


					locusTags.append(locus_tag)

					if locus_tag in filteredCDS[acc]:

						# Get the gene
						if record.features[feature].qualifiers.get('gene'):
							gene = record.features[feature].qualifiers.get('gene')
						else:
							gene = ''

						# Get the location
						location = str(record.features[feature].location)

						# Check to see whether there are multiple parts
						joincheck = re.search('join', location)

						if joincheck:

							joinCDSstart = 0
							joinCDSend = 0

							geneSeq = ''

							# Locate the splits
							region = location[location.find("{")+1:location.find("}")]


							splits = re.sub(', ', ',', region)
							splits = re.split(',', splits)

							cdsstrand = record.features[feature].strand


							intronNum = 0
							loc = ''

							# For each intron
							for i in range(0,len(splits)):

								intronNum += 1

								strand = record.features[feature].strand
								# strand = splits[i][splits[i].find("(")+1:splits[i].find(")")]

								locations = re.findall('\d+', splits[i])

								cdsStart = int(locations[0])
								cdsEnd = int(locations[1])


								if joinCDSstart == 0:
									joinCDSstart = int(cdsStart)
								else:
									if int(cdsStart) < joinCDSstart:
										joinCDSstart = int(cdsStart)

								if joinCDSend == 0:
									joinCDSend = int(cdsEnd)
								else:
									if int(cdsEnd) > cdsEnd:
										joinCDSend = int(cdsEnd)

								if intronNum == len(splits):
									loc += "%s..%s" % (cdsStart+1, cdsEnd)
								else:
									loc += "%s..%s," % (cdsStart+1, cdsEnd)



								seq = raw_sequence[cdsStart:cdsEnd]

								if strand == -1:
									strandType = 1
									seq = seq.reverse_complement()
									geneSeq += seq
								else:
									strandType = 0
									geneSeq += seq




							cdsStart = joinCDSstart
							cdsEnd = joinCDSend




						else:

							cdsStart = record.features[feature].location.nofuzzy_start
							cdsEnd = record.features[feature].location.nofuzzy_end

							loc = "%s..%s" % (cdsStart+1,cdsEnd)

							seq = raw_sequence[cdsStart:cdsEnd]

							strand = record.features[feature].strand


							if strand == -1:
								strandType = 1
								geneSeq = reverse_complement(seq)
							else:
								strandType = 0
								geneSeq = seq


						cds[acc][locus_tag] = geneSeq


	return cds



def compare_cds(leaders, filteredCDS, leaderFilePath):


	leaderGene = {}
	noLeaderGene = {}

	for acc in filteredCDS:



		leaderGene[acc] = [0,0]
		noLeaderGene[acc] = [0,0]

		for locus in filteredCDS[acc]:

			if locus in leaders[acc]:

				seq = filteredCDS[acc][locus]

				if seq[3] == "A":
					leaderGene[acc][0] += 1
				else:
					leaderGene[acc][1] += 1

			else:

				seq = filteredCDS[acc][locus]

				if seq[3] == "A":
					noLeaderGene[acc][0] += 1
				else:
					noLeaderGene[acc][1] += 1



	leaderFile = open(leaderFilePath, "w")

	leaderFile.write("acc,leader_a,leader_not_a,no_leader_a,no_leader_not_a\n")

	for acc in filteredCDS:



		leaderA = leaderGene[acc][0]
		leaderNotA = leaderGene[acc][1]
		noLeaderA = noLeaderGene[acc][0]
		noLeaderNotA = noLeaderGene[acc][1]

		totalLeader = float(leaderA + leaderNotA)
		totalNoLeader = float(noLeaderA + noLeaderNotA)

		fileLine = "%s,%s,%s,%s,%s\n" % (acc, leaderA/totalLeader, leaderNotA/totalLeader, noLeaderA/totalNoLeader, noLeaderNotA/totalNoLeader)
		leaderFile.write(fileLine)


#################################

def main():

	table4_genomes = get_table4_genomes(table4_genomes_path)

	# Retrieve the leader genes from the file
	leaders = get_leaders(filteredLeadersFile, table4_genomes)

	# Get the filtered genes
	filteredCDSloc = getFilteredGenes()

	# Get the distances from the CDS
	filteredCDS = get_filtered_cds(leaders, filteredCDSloc)



	compare_cds(leaders, filteredCDS, leaderFilePath)

#################################

if __name__ == "__main__":
	main()
