#!/usr/bin/python

# Script number: 			10.1
# File: 					1 of 4
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Identify potential leader genes residing upstream from the CDS
#							Trained using the methods from Korolev et al (2016) (https://biologydirect.biomedcentral.com/articles/10.1186/s13062-016-0123-8)
# Output file(s):			_leaderGenes.csv



import sys, os, re, shutil, Bio, csv
from sys import argv

from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate




##########################
# VARIABLES
##########################



# Set the directory with the raw good files
raw_genome_file_dir = 'bacterial_genomes/good_genomes/'
leader_directory = "outputs/leader_genes/"

testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["CP000925"]



##########################
# FUNCTIONS
##########################

# Setup new directories
def setupDirectory(directory):

	if not os.path.exists(directory):
		os.makedirs(directory)
		print "Created the new directory: %s" % directory

# Create a list of all the good files
def getGenomes(direc):

	files = []

	# For each file in the downloaded genomes
	for eachfile in os.listdir(direc):

		# if eachfile .startswith(testFile[0]):

		# If the file ends with the known .embl format
		if eachfile.endswith(".embl"):

			# Append the file to the files list
			files.append(eachfile)

	return files, len(files)


def newFile():



	# Output file with genes passing the filter
	leadersFilePath = leader_directory + "_leaderGenes.csv"
	leadersFile = open(leadersFilePath, "w")



	return leadersFile


def fileHead(file):

	headLine = "acc,locus_tag,gene_coords,leader_gene_coords\n"
	file.write(headLine)


def parseFile(acc,genome,fileNumber,numFiles, leadersFile):



	# Set the path to the good genome
	filePath = raw_genome_file_dir + genome

	# Get the acc number of the genome
	# acc = re.findall('.+(?=.embl)', genome)[0]

	# acc = "NC_002695"

	# # Set the path to the new file which is to be written
	# newPath = leader_directory + acc + ".txt"

	# # Open the new file
	# cdsFile = open(newPath, "w")

	# Read the file using BioPython
	# record = SeqIO.read(filePath, "gb")

	# print genome
	record = SeqIO.read(filePath, "embl")

	locusTags = []

	# Print miscellaneous information to the terminal
	print "-" * 10
	print record.name
	print record.description
	print "%s of %s" % (fileNumber, numFiles)






	for seq_record in SeqIO.parse(filePath, "embl"):

		raw_sequence = seq_record.seq



	cds = {}





	cdsNum = 0

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


			locusTags.append(locus_tag)

			# Get the gene
			if record.features[feature].qualifiers.get('gene'):
				gene = record.features[feature].qualifiers.get('gene')
			else:
				gene = ''

			# Get the location
			location = str(record.features[feature].location)

			# Check to see whether there are multiple exons
			joincheck = re.search('join', location)

			if joincheck:

				geneSeq = ''

				# Locate the exons
				region = location[location.find("{")+1:location.find("}")]


				exons = re.sub(', ', ',', region)
				exons = re.split(',', exons)

				intronNum = 0
				loc = ''

				# For each intron
				for i in range(0,len(exons)):

					intronNum += 1

					strand = exons[i][exons[i].find("(")+1:exons[i].find(")")]

					locations = re.findall('\d+', exons[i])

					cdsStart = int(locations[0])
					cdsEnd = int(locations[1])

					if intronNum == len(exons):
						loc += "%s..%s" % (cdsStart+1, cdsEnd)
					else:
						loc += "%s..%s," % (cdsStart+1, cdsEnd)



					seq = raw_sequence[cdsStart:cdsEnd]

					if strand == "-":
						strandType = 1
						seq = reverse_complement(seq)
						geneSeq += seq
					else:
						strandType = 0
						geneSeq += seq


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


			cds[cdsNum] = [strand, cdsStart,cdsEnd, locus_tag, geneSeq]



	numCDS = len(cds)
	leaders = {}

	identified = {}



	for gene in cds:


		strand = cds[gene][0]


		# If the gene is on the complement strand
		if strand == -1:


			if gene == numCDS:



				intergene = raw_sequence[cds[gene][2]:]
				geneSeq = raw_sequence[cds[gene][1]:cds[gene][2]]

				cdsStart = cds[gene][1]
				cdsEnd = cds[gene][2]

				intergeneStart = cds[gene][2]

				if len(intergene) > 1400:
					intergeneEnd = cdsEnd + 1400
				else:
					intergeneEnd = len(raw_sequence)



			else:

				intergene = raw_sequence[cds[gene][2]:cds[gene+1][1]]
				geneSeq = raw_sequence[cds[gene][1]:cds[gene][2]]

				intergene = intergene.reverse_complement()
				geneSeq = geneSeq.reverse_complement()

				# cdsStart is the 3' end of the CDS on the complement strand
				# cdsEnd is the 5' end of the CDS on the complement strand
				cdsStart = cds[gene][1]
				cdsEnd = cds[gene][2]

				# Intergene start is the 3' end of the intergenic region to the next CDS
				# Intergene end is the start of the next CDS
				intergeneStart = cds[gene][2]

				if len(intergene) > 1400:
					intergeneEnd = cdsEnd + 1400
				else:
					intergeneEnd = cds[gene+1][1]




			cdsregion = "%s..%s" % (cdsStart, cdsEnd)
			leaders[cdsregion] = {}

			leaders = find_orf(cdsregion, leaders, intergene, geneSeq, -1, translation_table)
			longestLeader = longest_orf(cdsregion, leaders, cdsStart, cdsEnd, intergeneStart, intergeneEnd, -1)

			# if cdsregion == "1694406..1694607":

			# 	print intergeneStart
			# 	print intergeneEnd
			# # 	print len(intergene)
			# 	# print intergene

			# 	print longestLeader



			# if gene == 15:
				# print "%s..%s" % (cdsStart, cdsEnd)
				# print "%s..%s" % (intergeneStart, intergeneEnd)

				# print cdsregion





				# print longestLeader[0], longestLeader[1]

				# Add 1 to the start to account for python index
				# print "%s..%s" % (intergeneEnd - longestLeader[0] - len(longestLeader[1]) + 1, intergeneEnd - longestLeader[0])

			# if cdsregion == "1549832..1551092":
			# 	print len(intergene)
			# 	print cdsregion, intergeneStart, intergeneEnd

			if len(longestLeader[1]) > 0:

				if len(intergene) > 1400:

					leaderLine = "%s,%s,complement(%s..%s),complement(%s..%s)\n" % (acc, cds[gene][3], cdsStart + 1, cdsEnd, intergeneStart + 1400 - longestLeader[0] - len(longestLeader[1]) + 1, intergeneStart + 1400 - longestLeader[0])
					leadersFile.write(leaderLine)
				else:
					leaderLine = "%s,%s,complement(%s..%s),complement(%s..%s)\n" % (acc, cds[gene][3], cdsStart + 1, cdsEnd, intergeneEnd - longestLeader[0] - len(longestLeader[1]) + 1, intergeneEnd - longestLeader[0])
					leadersFile.write(leaderLine)

		# If the gene is on the leading strand
		else:

			if gene == 1:

				intergene = raw_sequence[:cds[gene][1]]
				geneSeq = raw_sequence[cds[gene][1]:cds[gene][2]]

				cdsStart = cds[gene][1]
				cdsEnd = cds[gene][2]



				if len(intergene) > 1400:
					intergeneStart = cds[gene][1] - 1400
				else:
					intergeneStart = 0

				intergeneEnd = cds[gene][1]

			else:

				intergene = raw_sequence[cds[gene-1][2]:cds[gene][1]]
				geneSeq = raw_sequence[cds[gene][1]:cds[gene][2]]

				cdsStart = cds[gene][1]
				cdsEnd = cds[gene][2]

				if len(intergene) > 1400:
					intergeneStart = cds[gene][1] - 1400
				else:
					intergeneStart = cds[gene-1][2]

				intergeneEnd = cds[gene][1]



			cdsregion = "%s..%s" % (cdsStart, cdsEnd)
			leaders[cdsregion] = {}

			leaders = find_orf(cdsregion, leaders, intergene, geneSeq, 1, translation_table)

			longestLeader = longest_orf(cdsregion, leaders, cdsStart, cdsEnd, intergeneStart, intergeneEnd, 1)


			# if cdsregion =="5250..5547":
			# 	print longestLeader[1]
			# 	print longestLeader[1].translate()
			# if cdsregion == "2346841..2348446":

			# 	print intergene
			# 	print len(intergene)
			# 	print longestLeader[0]
			# 	print intergeneStart + longestLeader[0]
			# 	print longestLeader[1]




			if len(longestLeader[1]) > 0:

				# Add 1 to start to account for python indexing of the true sequence location
				# cds start, cds end, leader start, leader end (intergene start site + leader start site + length leader)
				leaderLine = "%s,%s,%s..%s,%s..%s\n" % (acc, cds[gene][3], cdsStart+1, cdsEnd, intergeneStart + longestLeader[0] + 1, intergeneStart + longestLeader[0] +  len(longestLeader[1]))
				leadersFile.write(leaderLine)




# Determine the longest ORF of potential leader genes
def longest_orf(cdsregion, leaders, cdsStart, cdsEnd, intergeneStart, intergeneEnd, strand):
	# Set the currest leader gene to empty
	# [leader start, leader seq, strand, cds start, cds end]
	longestLeader = [0, ""]

	# For each location of a potential leader gene
	for start in leaders[cdsregion]:

		# If there are no leader genes
		if longestLeader[1] == "":

			# Set the first possible leader gene
			longestLeader[0] = start
			longestLeader[1] = leaders[cdsregion][start]
		else:

			# If the ORF of the next potential leader gene is longer than the current longest, replace
			if len(leaders[cdsregion][start]) > len(longestLeader[1]):
				longestLeader[0] = start
				longestLeader[1] = leaders[cdsregion][start]

			# If the ORFs are the same length, choose the one closest to the CDS
			elif len(leaders[cdsregion][start]) == len(longestLeader[1]):

				# print "%s..%s - %snt " % (cdsStart, cdsEnd, cdsEnd - cdsStart)
				# print len(leaders[cdsregion][start]), len(longestLeader[1])

				if strand == 1:
					# Distance to the cds of the current longest ORF (cds start - leader gene end)
					currentLeaderDist = cdsStart - (longestLeader[0] + intergeneStart + len(longestLeader[1]))

					# Distance to the cds of the potential longest ORF (cds start - leader gene end)
					newLeaderDist = cdsStart - (intergeneStart + start + len(leaders[cdsregion][start]))

					# print "new: %s..%s - %s" % (intergeneStart + start, intergeneStart + start + len(leaders[cdsregion][start]), newLeaderDist)
					# print "current: %s..%s - %s" % (intergeneStart + longestLeader[0], intergeneStart + longestLeader[0] + len(longestLeader[1]), currentLeaderDist)

					if newLeaderDist < currentLeaderDist:
						longestLeader[0] = start
						longestLeader[1] = leaders[cdsregion][start]
				else:




					currentLeaderDist = (intergeneEnd - longestLeader[0] - len(longestLeader[1])) - cdsEnd
					newLeaderDist = (intergeneEnd - start - len(leaders[cdsregion][start])) - cdsEnd


					# print "%s..%s" % (cdsStart, cdsEnd)
					# print "Current: %s..%s - %snt" % (intergeneEnd - longestLeader[0] - len(longestLeader[1]) + 1, intergeneEnd - longestLeader[0], len(longestLeader[1]))
					# print "New: %s..%s - %snt" % (intergeneEnd - start - len(leaders[cdsregion][start]) + 1, intergeneEnd - start, len(leaders[cdsregion][start]))
					# print "-" * 20

					if newLeaderDist < currentLeaderDist:
						longestLeader[0] = start
						longestLeader[1] = leaders[cdsregion][start]

	return longestLeader




# Find all ORFs in the intergenic region
def find_orf(cdsregion, leaders, intergene, geneSeq, strand, trans_table):



	# Check the intergenic region is longer than 100 nucleotides long
	if len(intergene) >= 100:

		# if cdsregion == "2083736..2084318":
		# 		print intergene + "\n"

		# If the intergenic region is longer than 1400, only take
		# immediate upstream region of 1400 nucleotides
		if len(intergene) > 1400:

			# Select only the 1400 nucleotides upstream from the CDS
			intergene = intergene[-1400:]




		# if cdsregion == "4641570..4643259":
		# 	print intergene
		# 	print len(intergene)


		# Check that the CDS is between 200 and 100000 nucleotides long
		if len(geneSeq) > 200 and len(geneSeq) < 10000:

			# Start at each position in the intergenic region
			for i in range(0,len(intergene)):

				codon = intergene[i:i+3]

				# If the codon at the postion encodes a start codon
				if codon in ["ATG", "GTG"]:

					# Position in the intergene
					# print i
					orfseq = ""
					stop = False

					# For each subsequent codon
					for j in range(i, len(intergene),3):

						# Check if a stop codon has been encountered
						if stop == False:

							# If not, define the next codon
							codon = intergene[j:j+3]

							# If the codon is of the correct length
							if len(codon) % 3 == 0:

								# Add the codon to the sequence
								orfseq += codon

								# Check to see whether the codon is a stop codon

								if trans_table == "4":
									stop_codons = ["TAG, TAA"]
								else:
									stop_codons = ["TGA", "TAA", "TAG"]

								if codon in stop_codons:

									# If it is, set stop to true
									stop = True


					orflen = len(orfseq)

					# Check that the potential ORF is a mutile of 3
					if orflen % 3 == 0:

						# Ensure the potential ORF is longer than 5 codons
						if orflen > 18:

							# Ensure that the ORF ends with a regular stop codon
							if orfseq[-3:] in ["TGA", "TAA", "TAG"]:

								# print orfseq + "\n"
								leaders[cdsregion][i] = orfseq


	return leaders



#################################

def main():

	setupDirectory(leader_directory)
	leadersFile = newFile()
	fileHead(leadersFile)

	# Create the list of good files and number
	rawFiles, numFiles = getGenomes(raw_genome_file_dir)

	# rawFiles = ["bacterial_genomes/codons/ecoli_nc002695.gb"]

	fileNumber = 0




	for genome in rawFiles:



		acc = genome[:-5]

		#if acc in testFiles:

		fileNumber += 1
		parseFile(acc,genome,fileNumber,numFiles, leadersFile)



	leadersFile.close()

#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
