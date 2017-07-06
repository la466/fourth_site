#!/usr/bin/python

# Script number: 			14.4
# File: 					4 of 7
# Prerequisite script(s):	14.1, 14.2, 14.3
# Prerequisite file(s):		Raw embl genome files in good_genomes
# Description: 				Parse each of the good genomes and extract coding sequences (CDSs)



import sys, os, re, shutil, Bio, csv
from sys import argv

from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate

filename = argv



##########################
# VARIABLES
##########################


# Create the new directory to contain the genomes
newDir = "archaea_genomes/cds/"

# # Set the directory with the raw good files
goodDir = 'archaea_genomes/good_genomes/'


testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]
testFile = ["AE005174"]


##########################
# FUNCTIONS
##########################



# Set up the new directories for the good files and the bad files
def setupDirectories(direc):

	# If the directory already exists, delete and make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)


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




def parseFile(genome,fileNumber,numFiles):

	

	# Set the path to the good genome
	filePath = goodDir + genome

	# Get the acc number of the genome
	acc = re.findall('.+(?=.embl)', genome)[0]

	# Set the path to the new file which is to be written
	newPath = newDir + acc + ".txt"

	# Open the new file
	cdsFile = open(newPath, "w")
	
	# Read the file using BioPython
	record = SeqIO.read(filePath, "embl")

	locusTags = []

	# Print miscellaneous information to the terminal
	print "-" * 10
	print record.name
	print record.description
	print "%s of %s" % (fileNumber, numFiles)
	
	# Get the genus of the file
	with open(filePath) as myfile:

		DE = 0
		OC = 0

		for line in myfile:

			# IF the number of lines read starting with 'DE' is 0
			if DE == 0:

				# If the read lines starts with DE
				if line.startswith("DE"):

					# Increase the count of 'DE' so no more 'DE' lines will be read
					DE += 1

					#Extract the information from the line starting 'DE'
					reg = re.findall('(?<=DE\s\s\s).+', line)[0]

					# If the genus is contains 'UNVERIFIED'
					if 'UNVERIFIED' in reg:
						
						# Find the genus
						genus = re.findall('^UNVERIFIED:\s[a-zA-z0-9]+(?=\s)', reg)[0]

					else:

						# Find the genus
						genus = re.findall('([a-zA-Z0-9+]+)\s?', reg)[0]

			if line.startswith("OC"):
				
				OC += 1

				if OC == 1:
					reg = re.findall('(?<=OC\s\s\s).+', line)[0]
					splits = reg.split(';')
					phylum = splits[1]
					pattern = re.compile(r'\s+')
					phylum = re.sub(pattern, '', phylum)
			
	
	cdsNum = 0

	# Number of features for the file
	for feature in range(1,len(record.features)):

		# If the feature type is CDS
		if record.features[feature].type == "CDS":

			cdsNum += 1

			# Get the locus tag
			locus_tag = record.features[feature].qualifiers.get('locus_tag')
			translation_table = record.features[feature].qualifiers.get('transl_table')

			# If there is no locus tag, assign an undisclosed tag
			if not locus_tag:
				locus_tag = "no_tag_%s" % cdsNum
			else:
				locus_tag = locus_tag[0]

			
			samelocus = 0

			# Check to see if the locus is repeated
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

			# Check to see whether there are multiple cds_parts
			joincheck = re.search('join', location)


			# Check to see whether the cds has joins
			if joincheck:

				geneSeq = ''				

				# Locate the sequence
				region = location[location.find("{")+1:location.find("}")]

				
				cds_parts = re.sub(', ', ',', region)
				cds_parts = re.split(',', cds_parts)

				cds_part_num = 0
				loc = ''
				
				# For each sequence
				for i in range(0,len(cds_parts)):

					cds_part_num += 1

					# Get the strand the cds is on
					strand = cds_parts[i][cds_parts[i].find("(")+1:cds_parts[i].find(")")]				


					locations = re.findall('\d+', cds_parts[i])

					cdsStart = int(locations[0])
					cdsEnd = int(locations[1])

					if cds_part_num == len(cds_parts):
						loc += "%s..%s" % (cdsStart+1, cdsEnd)
					else:
						loc += "%s..%s," % (cdsStart+1, cdsEnd)



					seq = record.seq[cdsStart:cdsEnd]
					
					if strand == "-":
						strandType = 1
						seq = reverse_complement(seq)
						geneSeq += seq
					else:
						strandType = 0
						geneSeq += seq

				
			else:

				# Locate the sequence
				cdsStart = record.features[feature].location.nofuzzy_start
				cdsEnd = record.features[feature].location.nofuzzy_end

				loc = "%s..%s" % (cdsStart+1,cdsEnd)
				
				seq = record.seq[cdsStart:cdsEnd]

				strand = record.features[feature].strand				
				
				
				if strand == -1:
					strandType = 1
					geneSeq = reverse_complement(seq)
				else:
					strandType = 0
					geneSeq = seq


		
			# Write to file
			fileDesc = ">%s;phylum=%s;acc=%s;locus_tag=%s;gene=%s;loc=%s;complement=%s;trans_table=%s;gene_type=cds\n" % (genus,phylum,acc,locus_tag,gene,loc,strandType,translation_table[0])
			fileSeq = "%s\n\n" % geneSeq
			
		
			cdsFile.write(fileDesc+fileSeq)

	cdsFile.close()


		
#################################

def main():

	# Set up the new directory to print the genome extractions to
	setupDirectories(newDir)

	# Create the list of good files and number
	goodFiles, numFiles = getGenomes(goodDir)

	fileNumber = 0


	for genome in goodFiles:
		fileNumber += 1
		parseFile(genome,fileNumber,numFiles)


#################################

if __name__ == "__main__":
	main()

	



