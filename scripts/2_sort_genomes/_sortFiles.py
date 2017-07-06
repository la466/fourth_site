#!/usr/bin/python

# Script number:				2.1
# File:							1 of 1
# Prerequisite script(s):		1.1, 1.2
# Prerequisite file(s):		 	raw embl genome files (bacterial genome directory)
# Description: 					Sort the genomes into good and not needed genomes:
# 									1. Only include 1 genome per genus.
#									2. Genome must be over 500000 base pairs
#								Copy the good genomes to the new directory.

import sys, os, re, shutil

from sys import argv
from Bio import SeqIO


##########################
# VARIABLES
##########################


# Set the directory containing the downloaded raw genome files
genome_directory = 'bacterial_genomes/bacteria_raw_embl/'


# Set the new directories where the sorted files will be copied
newGood = 'bacterial_genomes/good_genomes/'
newBad = 'bacterial_genomes/bad_genomes/'
t4_genomes = 'bacterial_genomes/t4_genomes/'

good_genome_directory_output = 'good_genomes'




##########################
# FUNCTIONS
##########################

# Set up the new directories for the good files and the bad files
def setupDirectories(direc):

	# If the directory already exists, delete and make new directory
	if os.path.exists(direc):
		shutil.rmtree(direc)
	os.makedirs(direc)


# Create a list of all the files from the web crawler
def getGenomes(direc):

	files = []

	# For each file in the downloaded genomes
	for eachfile in os.listdir(direc):

		# If the file ends with the known .embl format
		if eachfile.endswith(".embl"):

			# Append the file to the files list
			files.append(eachfile)

	return files, len(files)



# Sort each of the genomes and copy to respective directory
def filterGenomes(genome,goodFiles,badFiles, output_file):

	# Set the path for the genome to be analysed
	genomePath = genome_directory +"/"+ genome

	# Open the file
	with open(genomePath) as myfile:

		table_find = False
		table_4 = False

		# Set the count of the file annotations to 0
		ID = 0
		DE = 0

		# Clean up the acc number of the genome
		name = re.findall('.+(?=.embl)', genome)[0]

		# For each line in the file being analysed
		for line in myfile:

			# If the file annotation is 'ID' and there have been no
			# previous occurraces of "ID"
			if line.startswith("ID") and ID == 0:

				# Extract the number of base pairs
				regex = re.findall('\d+\sBP', line)[0]

				# Extract the integer value of base pairs
				BP = int(re.findall('\d+', regex)[0])

				# Increase the count of 'ID' so no more 'ID' lines will be read
				ID += 1


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

			if not table_find:
				if line.startswith('FT'):
					reg = re.findall(r'(?<=FT)\s+\/transl_table=(\d+)', line)
					if len(reg) > 0:
						table = int(reg[0])
						table_find == True
						if table == 4:
							table_4 = True



		# If the genome length is >500000BP
		if BP > 500000:

			# Copy t4 genomes to directory
			if table_4:
				shutil.copy(genomePath, t4_genomes)

			# If the genus is not already used
			if genus not in goodFiles:

				# Add the genus to the goodFiles dictionary
				goodFiles[genus] = genome

				# Copy the file to the new good files directory
				shutil.copy(genomePath,newGood)
				shutil.copy(genomePath,good_genome_directory_output)

				print "Copying %s (%s) to good files (%s BP)" % (genus,goodFiles[genus],BP)
				output_file.write('%s,%s,%s\n' % (genome[:-5], genus, BP))



			# If the genus has already been used
			else:

				# Add the genome to the badFiles list
				badFiles.append(genome)

				# Copy the genome to the bad files directory
				shutil.copy(genomePath,newBad)




		# If the genome is <500000BP
		else:

			# Add the genome to the bad files list
			badFiles.append(genome)

			# Copy the genome to the bad files directory
			shutil.copy(genomePath,newBad)



def setupNiceDirectories(direc):

	# If the directory doesn't already exists, make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)


#################################

def main():

	# Set the number of files, good files and bad files to 0
	fileNumber = 0
	goodNum = 0
	badNum = 0

	# Create the dictionary and list to hold th good and bad files
	goodFiles = {}
	badFiles = []

	# Set up the new directories to copy the files to
	setupDirectories(newGood)
	setupDirectories(newBad)
	setupDirectories(t4_genomes)
	setupDirectories(good_genome_directory_output)


	setupNiceDirectories('outputs/')
	setupNiceDirectories('outputs/genome_info/')

	output_file = open('outputs/genome_info/genome_list.csv', 'w')
	output_file.write('accession,genus,base_pairs\n')

	# Return the list of downloaded genomes and the number
	genomes, numGenomes = getGenomes(genome_directory)

	# For each downloaded genome
	for genome in genomes:

		fileNumber +=1

		if fileNumber:
			print '-' * 20
			print genome
			print "%s of %s files" % (fileNumber,numGenomes)

			# Filter the genome on size and whether the genus has
			# already been used
			filterGenomes(genome,goodFiles,badFiles, output_file)

	output_file.close()

#################################


if __name__ == "__main__":
	main()
