#!/usr/bin/python

# Script number: 			10.3
# File: 					3 of 4
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 10.1, 10.3
# Prerequisite file(s):		_leader_genes_filtered.csv, _filteredGenes.csv, table4genomes.txt
# Description:				Analayse the data on leader genes
# 							1. Determine the total number of leader genes for each distance
# 							2. Calculate the propotion of CDS with +4 A given the CDS has a leader
# 							3. Comapre the GC3 content with the mean leader gene distance for each genome
# Output file(s):			_leader_distances.csv, _leader_fourth_site.csv, leader_distances_GC.csv, genomeLeaderDistances.csv, _genomeAPropLeaderDistances.csv



import re

from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate

import re


############
# Variables
############

output_dir = "outputs/leader_genes/"
filteredLeadersFile = output_dir + "_leader_genes_filtered.csv"
distanceFilePath = output_dir + "_leader_distances.csv"
emblPath = "bacterial_genomes/bacteria_raw_embl/"
leaderFourthPath = output_dir + "_leader_fourth_site.csv"
leaderDistanceGCPath = output_dir + "_leader_distances_GC.csv"
fourthSites = "outputs/ratio_testing/site_4/_site_4_ratios.csv"
table4_genomes_path = "outputs/gene_filtering/table4genomes.txt"

genomeLeaderDistanceAProp = output_dir + "_genomeLeaderDistances.csv"
genomeLeaderPropA_path = "outputs/leader_genes/_genomeAPropLeaderDistances.csv"





testFile = ["AE005174", "AE000511"]

############
# Functions
############




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


def leader_distances(leaders, distanceFilePath, gc3s, gcs, leaderDistanceGCPath, genomeLeaderDistanceAProp):

	distances = {}
	fourths = {}
	distancesGC = {}

	genomeDistances = {}


	genomeNo = 0
	numGenomes = len(leaders)



	for acc in leaders:
		# if acc in testFile:
		genomeNo += 1

			print acc
			print "%s of %s" % (genomeNo, numGenomes)
			print "-" * 20


			emblFile = emblPath + acc + ".embl"

			for seq_record in SeqIO.parse(emblFile, "embl"):
				raw_sequence = seq_record.seq





			distances[acc] = {}
			fourths[acc] = {}
			genomeDistances[acc] = {}


			# # Set up the dictionary containing distances
			# # for A content of CDS at each distance
			# for i in range(1,21):
			# 	genomeDistances[acc][i] = []


			if acc in gc3s:
				distancesGC[acc] = []
				gc3 = gc3s[acc]
				gc = gcs[acc]
				distancesGC[acc] = [gc3, [], gc]




			# For each CDS with a leader gene
			for locus in leaders[acc]:

				cds = leaders[acc][locus][0]
				leader = leaders[acc][locus][1]

				if "complement" in cds:

					# Minus 1 for python indexing
					cdsStart = int(re.findall("(?<=complement\()(\d+)(?=..)", cds)[0]) - 1
					cdsEnd = int(re.findall("(?<=\..)(\d+)(?=\)$)", cds)[0])

					leaderStart = int(re.findall("(?<=complement\()(\d+)(?=..)", leader)[0])
					leaderEnd = int(re.findall("(?<=\..)(\d+)(?=\)$)", leader)[0])

					leaderDistance = leaderStart - cdsEnd


					cdsSeq = raw_sequence[cdsStart:cdsEnd].reverse_complement()
					pos4 = cdsSeq[3]




				else:

					# Minus 1 for python indexing
					cdsStart = int(re.findall("^(\d+)(?=\..)", cds)[0]) - 1
					cdsEnd = int(re.findall("(?<=\..)(\d+)$", cds)[0])

					leaderStart = int(re.findall("^(\d+)(?=\..)", leader)[0])
					leaderEnd = int(re.findall("(?<=\..)(\d+)$", leader)[0])

					leaderDistance = (cdsStart+1) - leaderEnd


					cdsSeq = raw_sequence[cdsStart:cdsEnd]
					pos4 = cdsSeq[3]



				if leaderDistance in distances[acc]:
					distances[acc][leaderDistance].append(locus)
				else:
					distances[acc][leaderDistance] = [locus]

				# Get the fourth sites cds if there is a leader gene
				if leaderDistance in fourths[acc]:
					fourths[acc][leaderDistance].append(pos4)
				else:
					fourths[acc][leaderDistance] = [pos4]


				if acc in gc3s:
					distancesGC[acc][1].append(leaderDistance)


				if leaderDistance in genomeDistances[acc]:
					genomeDistances[acc][leaderDistance].append(pos4)
				else:
					genomeDistances[acc][leaderDistance] = [pos4]




	# Calculate the proprition of CDS with +4 with a distance under
	# the parameter i, given that a CDS has a leader gene
	genomeLeaderProps = {}

	for acc in genomeDistances:

		genomeLeaderProps[acc] = {}

		maxDistance = 0
		for leaderLength in genomeDistances[acc]:
			if leaderLength > maxDistance:
				maxDistance = leaderLength

		# print acc


		for i in range(1,21):

			genomeLeaderProps[acc][i] = [0,0]

			# Get the prop A for distances < i
			fourthAlower = 0
			totallower = 0

			for j in range(1,i+1):

				if j in genomeDistances[acc]:

					subsetA = genomeDistances[acc][j].count("A")
					lensubset = len(genomeDistances[acc][j])

					fourthAlower += subsetA
					totallower += lensubset

			# print "%s: %s, %s\n" % (i, fourthAlower, totallower)

			if fourthAlower == 0 or totallower == 0:
				genomeLeaderProps[acc][i][0] = 0
			else:
				fourthAlowerProp = fourthAlower / float(totallower)
				genomeLeaderProps[acc][i][0] = fourthAlowerProp



			# Get the prop A for distances < i
			fourthAupper = 0
			totalupper = 0

			for j in range(i+1,maxDistance+1):
				if j in genomeDistances[acc]:

					subsetAupper = genomeDistances[acc][j].count("A")
					lensubsetupper = len(genomeDistances[acc][j])

					fourthAupper += subsetAupper
					totalupper += lensubsetupper

			# print "%s: %s, %s\n" % (i, fourthAupper, totalupper)

			if fourthAupper == 0 or totalupper == 0:
				genomeLeaderProps[acc][i][1] = 0
			else:
				fourthAupperProp = fourthAupper / float(totalupper)
				genomeLeaderProps[acc][i][1] = fourthAupperProp


	genomeLeader = open(genomeLeaderDistanceAProp, "w")
	fileLine = "acc,total_leaders"
	for i in range(1,21):
		fileLine += ",under_%s_a_prop,over_%s_a_prop" % (i,i)
	fileLine += "\n"
	genomeLeader.write(fileLine)

	for acc in genomeLeaderProps:

		fileLine = "%s,%s" % (acc,len(leaders[acc]))

		for distance in genomeLeaderProps[acc]:


			fileLine += ",%s,%s" % (genomeLeaderProps[acc][distance][0], genomeLeaderProps[acc][distance][1])
			# print acc, distance, genomeLeaderProps[acc][distance][0], genomeLeaderProps[acc][distance][1]
		fileLine += "\n"

		genomeLeader.write(fileLine)
	genomeLeader.close()



	#Get the proportion of A content at each distances for each genome
	genomeLeaderPropADistances = {}
	for acc in genomeDistances:

		genomeLeaderPropADistances[acc] = {}


		for i in range(1,52):

			genomeLeaderPropADistances[acc][i] = 0

			# Get the prop A for distances < i
			fourthACount = 0
			fourthTotal = 0


			if i in genomeDistances[acc]:

				fourthACount = genomeDistances[acc][i].count("A")
				fourthTotal = len(genomeDistances[acc][i])


			if fourthACount == 0 or fourthACount == 0:
				genomeLeaderPropADistances[acc][i] = 0
			else:
				genomeLeaderPropADistances[acc][i] = fourthACount / float(fourthTotal)



	genomeLeaderPropA = open(genomeLeaderPropA_path, "w")
	fileLine = "acc,total_leaders"
	for i in range(1,51):
		fileLine += ",%s" % (i)
	fileLine += "\n"
	genomeLeaderPropA.write(fileLine)

	for acc in genomeLeaderPropADistances:
		fileLine = "%s" % acc

		for i in range(1,52):
			fileLine += ",%s" % genomeLeaderPropADistances[acc][i]

		fileLine += "\n"

		genomeLeaderPropA.write(fileLine)

	genomeLeaderPropA.close()
















	# Calculate the nucleotide content of the fourth site
	# given a CDS has a leader gene
	maxDistance = 0
	meanDistances = {}

	for acc in distances:
		for distance in sorted(distances[acc]):
			if distance > maxDistance:
				maxDistance = distance


	allFourths = {}
	for i in range(1, maxDistance +1):
		allFourths[i] = []


	for acc in fourths:
		for distance in sorted(fourths[acc]):

			allFourths[distance].extend(fourths[acc][distance])


	leaderFourthFile = open(leaderFourthPath, "w")
	fileLine = "distance,a_content,c_content,t_content,g_content\n"
	leaderFourthFile.write(fileLine)

	for distance in allFourths:

		if len(allFourths[distance]) == 0:
			acontent, ccontent, tcontent, gcontent = [0,0,0,0]
		else:
			acontent = allFourths[distance].count("A") / float(len(allFourths[distance]))
			ccontent = allFourths[distance].count("C") / float(len(allFourths[distance]))
			tcontent = allFourths[distance].count("T") / float(len(allFourths[distance]))
			gcontent = allFourths[distance].count("G") / float(len(allFourths[distance]))

		fileLine = "%s,%s,%s,%s,%s\n" % (distance, acontent, ccontent, tcontent, gcontent)
		leaderFourthFile.write(fileLine)





	# Open the file and write the header
	leaderFile = open(distanceFilePath, "w")

	leaderFile.write('acc')
	for i in range(1, maxDistance + 1):
		headerRange = ",%s" % i
		leaderFile.write(headerRange)
	leaderFile.write("\n")



	for i in range(1, maxDistance + 1):
		meanDistances[i] = []


	# For each of the genomes
	for acc in distances:
		leaderFile.write(acc)

		# For each distance from the CDS
		for i in range(1, maxDistance + 1):

			# Add the number of leader genes that distance from the CDS
			if i in distances[acc]:
				fileLine = ",%s" % len(distances[acc][i])
				leaderFile.write(fileLine)

				meanDistances[i].append(len(distances[acc][i]))

			else:
				leaderFile.write(",0")





		leaderFile.write("\n")



	leaderFile.write("total")
	for distance in meanDistances:
		# print distance, sum(meanDistances[distance])
		fileLine = ",%s" % sum(meanDistances[distance])
		leaderFile.write(fileLine)
	leaderFile.close()





	# Create a file showing the mean leader gene distance and
	# GC3 content to see whether there is a correlation between the two
	distancesGCFile = open(leaderDistanceGCPath, "w")
	distancesGCFile.write("acc,gc3,gc,mean_leader_dist\n")

	for acc in distancesGC:
		meanDistance = sum(distancesGC[acc][1]) / float(len(distancesGC[acc][1]))
		fileLine = "%s,%s,%s,%s\n" % (acc,distancesGC[acc][0],distancesGC[acc][2], meanDistance)
		distancesGCFile.write(fileLine)

	distancesGCFile.close()



def get_GC(fourthSites):

	GC3s = {}
	GCs = {}
	rowNum = 0

	with open(fourthSites) as myfile:

		for row in myfile:

			rowNum += 1

			if rowNum > 1:

				genus = row.split(",")

				acc = genus[2]
				gc3 = genus[6]
				gc = genus[5]

				GC3s[acc] = gc3
				GCs[acc] = gc



	return GC3s, GCs



#################################

def main():

	# Get GC3 content for each genome
	gc3, gc = get_GC(fourthSites)

	table4_genomes = get_table4_genomes(table4_genomes_path)

	# Retrieve the leader genes from the file
	leaders = get_leaders(filteredLeadersFile, table4_genomes)


	# Get the distances from the CDS
	leader_distances(leaders, distanceFilePath, gc3, gc, leaderDistanceGCPath, genomeLeaderDistanceAProp)


#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
