#!/usr/bin/python

# Script number: 			13.1
# File: 					1 of 3
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Get the 3' tails for the genomes
# Output file(s):			antiSDsequences.txt


import sys, os, re, shutil, Bio, csv, random, subprocess, urllib, collections
from sys import argv

from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from collections import Counter





##########################
# VARIABLES
##########################


# Set the directory with the raw good files
rawFileDir = 'bacterial_genomes/bacteria_raw_embl/'
filtered_genes_file_path = "outputs/gene_filtering/_filteredGenes.csv"
t4genomesPath = "outputs/gene_filtering/table4genomes.txt"
antiSDfilePath = "outputs/sd_sequences/antiSDsequences.txt"
outputDir = "outputs/"
output_sd_dir = outputDir + "sd_sequences"

# testFiles = ["AE000657", "AE014295", "AE005174", "BA000019", "AL591688", "AP006840", "AP012205", "AE000512", "AE017221", "AE009952"]
testFiles = ["AE005174", "AE000657", "AE000511", "AE003852", "AE004091", "AE0064648", "AE008691", "AE003849", "AL111168"]
testFiles2 = ["AE000657", "AE014295", "AE005174", "AE017198", "BA000019", "AL591688", "AP006840", "AP012205", "AE000512", "AE017221", "AE009952", "AE006840"]
testFile = ["AE005174"]



##########################
# FUNCTIONS
##########################



def remove_directory(directory):
	if os.path.exists(directory):
		shutil.rmtree(directory)


def remove_make_directory(directory):

	if os.path.exists(directory):
		shutil.rmtree(directory)
	os.makedirs(directory)


# Set up the new directories for the good files and the bad files
def setupDirectory(direc):

	# If the directory already exists, delete and make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)



def get_T4_genomes(T4path):

	T4genomes = []

	with open(T4path, "U") as myfile:

		for line in myfile:

			line = line.strip("\n")
			splits = line.split(",")

			T4genomes.append(splits[1])

	return T4genomes



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



def getFilteredGenes(filtered_genes_file_path):


	filteredGenes = {}


	with open(filtered_genes_file_path, 'U') as myfile:
		read = csv.reader(myfile)



		rowNum = 0
		for row in read:
			rowNum += 1
			if rowNum > 1:


				acc = row[1]

				filteredGenes[acc] = []


				for i in range(2, len(row)):
					filteredGenes[acc].append(row[i])




	return filteredGenes



def get_sequence(feature, raw_sequence, record, locusTags, rRNAnum, trans_table):


	# Get the locus tag
	locus_tag = record.features[feature].qualifiers.get('locus_tag')

	# Set the locus tag
	if not locus_tag:
		locus_tag = "no_tag_%s" % rRNAnum
	else:
		locus_tag = locus_tag[0]


	# Check to see whether the locus has already been used
	if locus_tag in locusTags:
		locus_tag = "%s_%s" % (locus_tag, samelocus + 1)
		quit()
	else:
		locusTags.append(locus_tag)




	# Get the location of the gene
	# Note, the start site is -1 to the true location
	location = str(record.features[feature].location)


	# Check to see whether there are multiple parts
	joincheck = re.search('join', location)

	# If there is a join in the gene
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



	# Otherwise
	else:

		# Locate the gene start (note -1 from true start site)
		cdsStart = record.features[feature].location.nofuzzy_start

		# Locate the gene end
		cdsEnd = record.features[feature].location.nofuzzy_end

		# Write the location
		loc = "%s..%s" % (cdsStart+1,cdsEnd)

		# Get the raw sequence
		seq = raw_sequence[cdsStart:cdsEnd]



		strand = record.features[feature].strand


		if strand == -1:
			strandType = 1
			geneSeq = reverse_complement(seq)
		else:
			strandType = 0
			geneSeq = seq






	return geneSeq, cdsStart, cdsEnd, strand


def get_RNA_tails(record, raw_sequence):



	# Retrieve the 16S rRNA for each genomes
	print "...Extracting rRNAs"


	rRNAtails = {}
	locusTags = []
	rRNAnum = 0
	trans_table = False

	# Foreach feature in the genome
	for feature in range(1,len(record.features)):



		# Get the translation table for the genome
		translation_table = record.features[feature].qualifiers.get('transl_table')
		if translation_table and trans_table == False:
			trans_table = translation_table[0]

		# If the feature type is an rRNA
		if record.features[feature].type == "rRNA":

			# Increase the rRNA count
			rRNAnum += 1

			# Extract the product / note of the feature
			product = record.features[feature].qualifiers.get("product")
			note = record.features[feature].qualifiers.get("note")

			# Set the 16S rRNA identifier to false
			rRNA16S = False

			# Check to see whether the feature is a 16S rRNA, if so set the
			# 16S rRNA identitier to true
			if product:
				if product[0] == "16S ribosomal RNA" or product[0] == "ribosomal RNA-16S" or product[0] == "16S rRNA":
					rRNA16S = True
			elif note:
				if note[0] == "16S rRNA":
					rRNA16S = True

			# If the 16 rRNA identifier is true
			if rRNA16S == True:

				# Extract the 16S rRNA sequence
				geneSeq, cdsStart, cdsEnd, strand = get_sequence(feature, raw_sequence, record, locusTags, rRNAnum, trans_table)

				loc = geneSeq.rfind("GAT")

				rRNAend = geneSeq[loc:]
				# print rRNAend

				# Add the tail to the dictionary
				if str(rRNAend) in rRNAtails:
					rRNAtails[str(rRNAend)] += 1
				else:
					rRNAtails[str(rRNAend)] = 1


	# # Get the tail that occurs the most number of times
	# # Retrieve the 16S rRNA for each genomes
	print "...Calculating the most common tail"
	if rRNAtails:
		maxTail = max(rRNAtails, key=rRNAtails.get)
	else:
		maxTail = False

	return maxTail




# Get all the CDS feature numbers
def get_cds_feature_numbers(record, cdsLocations):

	# Foreach feature in the genome
	for feature in range(1,len(record.features)):

		# If the feature type is an rRNA
		if record.features[feature].type == "CDS":

			cdsLocations.append(feature)


	return cdsLocations







def get_upstream_seq(record, raw_sequence, feature, filteredGenes, acc):

	cdsStart = record.features[feature].location.nofuzzy_start
	cdsEnd = record.features[feature].location.nofuzzy_end
	cdsstrand = record.features[feature].strand
	locus_tag = record.features[feature].qualifiers.get('locus_tag')


	if locus_tag:


		if locus_tag[0] in filteredGenes[acc]:

			if cdsstrand == -1:
				upstreamSeq = reverse_complement(raw_sequence[cdsEnd:cdsEnd+20])
			else:
				upstreamSeq = raw_sequence[cdsStart-20:cdsStart]

		else:
			locus_tag[0] = ""
			upstreamSeq = ""

	return upstreamSeq, locus_tag[0]




def get_rRNA_tail(acc,genome,fileNumber,numFiles, filteredGenes, maxTailsShort, maxTailsLong, antiSDfile, accessions):



	# Set the path to the good genome
	filePath = rawFileDir + genome

	record = SeqIO.read(filePath, "embl")

	# Print miscellaneous information to the terminal
	print "-" * 10
	print record.name
	print record.description
	print "%s of %s" % (fileNumber, len(accessions))


	acc = record.name



	# Get the raw sequence for the genomes
	print "...Getting raw sequence"
	for seq_record in SeqIO.parse(filePath, "embl"):
		raw_sequence = seq_record.seq

	# Get CDS feature numbders
	print "...Getting CDS feature numbers"
	cdsLocations = []
	get_cds_feature_numbers(record, cdsLocations)

	# Write the CDS sequences to a fasta file
	# write_cds_to_fasta(record, filteredGenes, acc, cdsLocations, raw_sequence)




	# Get the sequence of the most common 16S rRNA tail
	print "...Getting the most common 16S rRNA tail"
	maxTail = get_RNA_tails(record, raw_sequence)

	if maxTail != False:
		if len(maxTail) <= 15 and len(maxTail) >= 8:
			maxTailsShort.append(acc)

			fileLine = ">%s\n" % acc
			fileLine += "%s\n" % maxTail
			antiSDfile.write(fileLine)

	else:
			maxTailsLong.append(acc)




def retrieve_rRNA_tails(filePath):

	if not os.path.isfile(filePath):
		print "No anti Shine Dalgarno sequences defined. Use '_sdSequences.py get' to retrieve"
		exit()


	antiSDsequences = {}

	with open(filePath, "U") as myfile:

		lines = myfile.read()

		# Locate each of the genomes
		genomes = re.findall('>.*\n.*(?=\n)', lines)

		for genome in genomes:

			# Locate the sequence
			acc = re.findall('(?<=>).*(?=\n)', genome)[0]
			antiSDseq = re.findall('\n(.*)', genome)[0]
			antiSDsequences[acc] = antiSDseq

	return antiSDsequences



def retrieve_upsstream_sequences(acc,fileNumber,filteredGenes):

	upstreamCDS = {}

	# Set the path to the good genome
	filePath = rawFileDir + acc + ".embl"

	record = SeqIO.read(filePath, "embl")

	# Print miscellaneous information to the terminal
	print "-" * 10
	print acc
	print record.description
	print "Genome %s" % (fileNumber)

	trans_table = False
	cdsCount = 0

	# Get the raw sequence for the genomes
	print "...Getting raw sequence"
	for seq_record in SeqIO.parse(filePath, "embl"):
		raw_sequence = seq_record.seq

	# Foreach feature in the genome
	for feature in range(1,len(record.features)):

		# Get the translation table for the genome
		translation_table = record.features[feature].qualifiers.get('transl_table')
		if translation_table and trans_table == False:
			trans_table = translation_table[0]

		# If the feature type is an rRNA
		if record.features[feature].type == "CDS":

			# Get the 20BP upstream
			upstreamSeq, locus_tag = get_upstream_seq(record, raw_sequence, feature, filteredGenes, acc)

			# If the cds is in the filtered genes
			if upstreamSeq and locus_tag:
				cdsCount +=1

			upstreamCDS[locus_tag] = upstreamSeq

	return upstreamCDS







def setup_free2bind(genome, antiSDsequences, genomeNum, filteredGenes,upstreamDirectory):


	genomeNum += 1


	# Set the path to the good genome
	filePath = rawFileDir + genome + ".embl"

	record = SeqIO.read(filePath, "embl")

	# Print miscellaneous information to the terminal
	print "-" * 10
	print genome
	print record.description
	print "Genome %s" % (genomeNum)

	trans_table = False
	cdsCount = 0

	# Get the raw sequence for the genomes
	print "...Getting raw sequence"
	for seq_record in SeqIO.parse(filePath, "embl"):
		raw_sequence = seq_record.seq

	antiSD = antiSDsequences[genome]



	upstreamGenomeDirectory = upstreamDirectory + genome + "/"
	setupDirectory(upstreamGenomeDirectory)

	# Write raw sequence to fasta file
	upstreamRawSeqPath = upstreamGenomeDirectory + genome + "_seq.fasta"
	rawSeqFile =  open(upstreamRawSeqPath, "w")
	fileLine = ">%s\n" % genome
	fileLine += "%s" % raw_sequence
	rawSeqFile.write(fileLine)
	rawSeqFile.close()

	# Get the locations or the upstream/downstream regions centered around the start codon

	# Setup file
	upstreamLocationsPath = upstreamGenomeDirectory + genome + "_index_data"
	upstreamLocationsFile = open(upstreamLocationsPath, "w")

	# Foreach feature in the genome
	for feature in range(1,len(record.features)):

		# Get the translation table for the genome
		translation_table = record.features[feature].qualifiers.get('transl_table')
		if translation_table and trans_table == False:
			trans_table = translation_table[0]

		# If the feature type is an rRNA
		if record.features[feature].type == "CDS":

			locus_tag = record.features[feature].qualifiers.get('locus_tag')
			if locus_tag:
				if locus_tag[0] in filteredGenes[genome]:

					# Get the strand the sequence is on
					strand = record.features[feature].strand

					# Locate the gene start (note -1 from true start site)
					cdsStart = record.features[feature].location.nofuzzy_start

					# Locate the gene end
					cdsEnd = record.features[feature].location.nofuzzy_end

					# Get the 60bp surrounding the start codon
					# plus then length of the antiSD to calcaulte to for the
					# last base in the 60bp region
					if strand == -1:
						loc = "complement(%s..%s)" % (cdsEnd-30-len(antiSD)+1,cdsEnd+30)
					else:
						loc = "%s..%s" % (cdsStart-30,cdsStart+30+len(antiSD)-1)

					fileLine = "%s # %s\n" % (loc, locus_tag[0])
					upstreamLocationsFile.write(fileLine)

	upstreamLocationsFile.close()


	print "...Calculcating dG values"
	print "...AntiSD: 3' - %s - 5'" % antiSD[::-1]

	currentDir = os.getcwd()
	relPathToLaunch =  currentDir + "/free2bind/"

	# Copy the genome files to the free2bind directory
	shutil.copy(upstreamLocationsPath, relPathToLaunch)
	shutil.copy(upstreamRawSeqPath, relPathToLaunch)

	# Change directory to the free2bind directory
	os.chdir(relPathToLaunch)

	# Get the input and output files
	fasta_file = genome + "_seq.fasta"
	index_file = genome + "_index_data"
	output_file = genome + "_dg_values.txt"
	output_dir = relPathToLaunch + "dg_values/"
	print output_dir

	setupDirectory(output_dir)

	# Run free2bind
	processCommand = "./free_scan.pl %s %s %s %s" % (antiSD[::-1], fasta_file, index_file, "dg_values/")
	subprocess.call(processCommand, shell=True)

	# Return to the previous directory
	os.chdir(currentDir)

	output_file_path = relPathToLaunch + output_file

	# Remove dG directory if it already exists
	outputFinalDest = upstreamGenomeDirectory + "dg_values/"
	remove_directory(outputFinalDest)

	# Return the dG values to the genome directory
	shutil.move(output_dir, upstreamGenomeDirectory)

	# Delete the files from the free2bind directory
	os.remove(relPathToLaunch + fasta_file)
	os.remove(relPathToLaunch + index_file)


def extract_dg_values(genome, antiSDsequences, genomeNum, upstreamDirectory, filteredGenes, dGvalues, mean_dG_values):



	dgFileDir = upstreamDirectory + genome + "/dg_values/"

	sdA_location = len(antiSDsequences[genome]) - antiSDsequences[genome].find("ACC")

	# If there a dG values
	if len(os.listdir(dgFileDir)) != 0:

		dGvalues[genome] = {}
		fileNo = 0

		print "...Extracting dG values for %s" % genome
		gene_sequences = get_gene_sequences_for_dg(genome, filteredGenes)


		# For each of the genes
		for file in os.listdir(dgFileDir):

			if '.txt' in file:

				dgFilePath = dgFileDir + file


				# Get the locus of the gene analysed
				locus = file.replace("_dg.txt", "")




				genome_dGvalues = {}

				with open(dgFilePath, 'U') as dGfile:

					fileNo += 1
					# if fileNo == 1:

					lines = dGfile.read()

					# Extract the dGs, removing the header information
					dGs = re.sub(r'(?m)^\#.*\n?', '', lines)

					# Split each individual dG value
					dGs = dGs.split('\n')

					# Remove the last 2 blank entries from \n
					dGs = dGs[:-2]

					# For each of the dG entries
					for dGentry in dGs:

						# Get the value
						dG_value = re.findall('^(.*)\\t', dGentry)
						dG_value = float(dG_value[0])

						# Get the binding position
						dG_binding_pos = re.findall('^.*\\t\#\s\w\s(\d*)$', dGentry)
						dG_binding_pos = int(dG_binding_pos[0])

						# print dG_binding_pos
						genome_dGvalues[dG_binding_pos] = dG_value

						# Append to the mean dG dictionary
						sdA_relative_location = dG_binding_pos-30 + sdA_location
						if sdA_relative_location in mean_dG_values:



							if locus in gene_sequences:
								mean_dG_values[dG_binding_pos-30 + sdA_location].append(dG_value)

					max_dG = [0,0]




					for pos in genome_dGvalues:

						if genome_dGvalues[pos] < max_dG[1]:
							max_dG = [pos-30, genome_dGvalues[pos]]


					# print locus
					# print gene_sequences[locus]
					# print max_dG
					if locus in gene_sequences:
						dGvalues[genome][locus] = [max_dG[0], max_dG[1], gene_sequences[locus]]



def get_gene_sequences_for_dg(genome, filteredGenes):


	gene_sequences = {}

	geneFilePath = codonDir + "genomes/" + genome + ".txt"

	with open(geneFilePath, "r") as geneFile:

		lines = geneFile.read()

		# Locate each of the strains
		genes = re.findall('>.*\n.*(?=\n)', lines)

		geneNum = 0

		# For each of the strains
		for gene in genes:

			geneNum += 1

			# Locate the strain identifier
			gene_locus = re.findall("(?<=locus_tag=)(.*)(?=;gene=)", gene)[0]
			# Locate the sequence
			seq = re.findall('\n(.*)', gene)[0]

			if gene_locus in filteredGenes[genome]:
				gene_sequences[gene_locus] = seq


	return gene_sequences


def analyse_dgs(dGvalues, upstreamDirectory, GC3s, Aratios, antiSDsequences):



	print "\n...Analysing dG values\n"


	genomeOutputFilePath = codonDir + "_sdSequenceAnalysis.csv"
	genomeOutputFile = open(genomeOutputFilePath, "w")
	genomeOutputFile.write('genome,gc3,a_ratio,num_genes,prop_sd,prop_no_sd,prop_sd_a,prop_no_sd_a,prop_weak_sd_a,prop_strong_sd_a,near_prop,far_prop\n')

	sdPosCount = {}
	for i in range(-30,30):
		sdPosCount[i] = {}

	for genome in dGvalues:

		if genome in os.listdir(upstreamDirectory):

			print genome
			bind = []
			noBind = []

			bindA = []
			noBindA = []

			weakBindA = []
			strongBindA = []


			sd_distances = []



			sdA_location = len(antiSDsequences[genome]) - antiSDsequences[genome].find("ACC")


			# Calcaulte the numbder of SD sequences at each position for each genome
			for pos in sdPosCount:
				sdPosCount[pos][genome] = 0


			for locus in sorted(dGvalues[genome]):
				distance = dGvalues[genome][locus][0]
				dgValue = dGvalues[genome][locus][1]


				if dgValue < -3.4535:

					if distance+sdA_location-1 < 0:
						sd_distances.append(distance+sdA_location-1)

					if distance+sdA_location in sdPosCount:

						# distance - length to 5' - ACC - 3' -1 to set A of ATG as 0
						sdPosCount[distance+sdA_location-1][genome] += 1


			# Get the most common upstream sd 5' - ACC - 3' A position

			most_common_dist = float(Most_Common(sd_distances))


			near = []
			far = []

			for locus in sorted(dGvalues[genome]):



				distance = dGvalues[genome][locus][0]
				dgValue = dGvalues[genome][locus][1]
				fourth = dGvalues[genome][locus][2][3]

				sd_distance = float(distance+sdA_location-1)


				if dgValue <= -3.4535 and sd_distance < 0:
					if sd_distance >= most_common_dist:
						near.append(fourth)
					else:
						far.append(fourth)



			if len(near) != 0:
				near_prop =  near.count("A") / float(len(near))
			else:
				near_prop = 0

			if len(far) != 0:
				far_prop = far.count("A") / float(len(far))
			else:
				far_prop = 0

			# Calculate whether an SD sequence exits, and the strength of SDs
			# and get the fourth base of the sequence
			for locus in sorted(dGvalues[genome]):

				distance = dGvalues[genome][locus][0]
				dgValue = dGvalues[genome][locus][1]
				fourth = dGvalues[genome][locus][2][3]

				sd_distance = float(distance+sdA_location-1)

				# Ensure the SD starts downstream
				if dgValue <= -3.4535 and sd_distance <= 30:
					bind.append(locus)
					bindA.append(fourth)

				else:
					noBind.append(locus)
					noBindA.append(fourth)

				if dgValue <= -3.4535 and dgValue >= -8.4 and sd_distance <= 30:
					weakBindA.append(fourth)
				elif dgValue <= -8.4 and distance <= 30:
					strongBindA.append(fourth)


			# Calculate proportions of A
			if len(bind) != 0 or len(noBind) != 0:

				prop_sd = len(bind) / float(len(bind) + len(noBind))
				prop_no_sd = len(noBind) / float(len(bind) + len(noBind))

				prop_a_sd = bindA.count("A") / float(len(bindA))
				prop_a_no_sd = noBindA.count("A") / float(len(noBindA))

				if len(weakBindA) != 0:
					prop_weak_sd_a = weakBindA.count("A") / float(len(weakBindA))
				else:
					weakBindA = 0

				if len(strongBindA) != 0:
					prop_strong_sd_a = strongBindA.count("A") / float(len(strongBindA))
				else:
					prop_strong_sd_a = 0

				fileLine = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (genome, GC3s[genome], Aratios[genome], len(dGvalues[genome]), prop_sd, prop_no_sd, prop_a_sd, prop_a_no_sd, prop_weak_sd_a, prop_strong_sd_a,near_prop,far_prop)
				genomeOutputFile.write(fileLine)



	# Output the number of genes at each distance
	outputFilePath = codonDir + "_sd_distances.csv"
	outputFile = open(outputFilePath, "w")

	genomes = []
	fileLine = "pos"
	for pos in sorted(sdPosCount):
		for genome in sorted(sdPosCount[pos]):
			if genome not in genomes:
				genomes.append(genome)

	for genome in genomes:
		fileLine += ",%s" % genome
	fileLine += "\n"
	outputFile.write(fileLine)

	for pos in sorted(sdPosCount):

		fileLine = "%s" % pos
		for genome in sorted(sdPosCount[pos]):
			fileLine += ",%s" % (sdPosCount[pos][genome])

		fileLine += "\n"
		outputFile.write(fileLine)

	outputFile.close()

def Most_Common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]


def get_accessions():

	accessions = []

	for file in os.listdir('genome_extractions/cds'):
		if file.endswith('.txt'):
			accessions.append(file.strip('.txt'))

	return (accessions)

#################################

def main():

	# Retrieve the T4 genomes
	T4genomes = get_T4_genomes(t4genomesPath)

	accessions = get_accessions()


	# Create the list of good files and number
	rawFiles, numFiles = getGenomes(rawFileDir)

	# rawFiles = ["bacterial_genomes/codons/ecoli_nc002695.gb"]

	fileNumber = 0

	# Get a list of the filtered genes
	filteredGenes = getFilteredGenes(filtered_genes_file_path)

	maxTailsLong = []
	maxTailsShort = []


	setupDirectory(output_sd_dir)


	#if "-get" in argv:
	# Determine anti SD sequences
	antiSDfile = open(antiSDfilePath, "w")
	for genome in rawFiles:



		acc = genome[:-5]
		# if acc in testFile:
		if acc not in T4genomes:
			if acc in accessions:


			# if acc in testFile:
				fileNumber += 1
				get_rRNA_tail(acc,genome,fileNumber,numFiles, filteredGenes, maxTailsShort, maxTailsLong, antiSDfile, accessions)

	antiSDfile.close()


	# # Retrieve the anti SD sequence
	# antiSDsequences = retrieve_rRNA_tails(antiSDfilePath)
	#
	# # Locate the upstream sequences for each genome with a
	# # suitable anti SD sequence
	# fileNumber = 0
	#
	#
	# upstreamDirectory = "bacterial_genomes/codons/antiSD_upstream_seq/"
	# setupDirectory(upstreamDirectory)
	#
	#
	# dGvalues = {}
	# mean_dG_values = {}
	# for i in range(-30,31):
	# 	mean_dG_values[i] = []
	#
	# for genome in sorted(antiSDsequences):
	#
	# 	# if genome in testFiles2:
	#
	#
	# 		if "-calc_dg" in argv:
	# 			genomeNum = 0
	# 			setup_free2bind(genome, antiSDsequences, genomeNum, filteredGenes, upstreamDirectory)
	#
	# 		if "-extract_dg" in argv:
	#
	# 			if genome in os.listdir(upstreamDirectory):
	# 				genomeNum = 0
	# 				extract_dg_values(genome, antiSDsequences, genomeNum, upstreamDirectory, filteredGenes, dGvalues, mean_dG_values)
	#
	# analyse_dgs(dGvalues, upstreamDirectory, GC3s, Aratios, antiSDsequences)
	#
	# for pos in sorted(mean_dG_values):
	# 	if len(mean_dG_values[pos]) != 0:
	# 		print pos, sum(mean_dG_values[pos]) / float(len(mean_dG_values[pos]))
	#
	#
	#



#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
