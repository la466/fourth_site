#!/usr/bin/python

# Script number:				14.2
# File:							2 of 7
# Prerequisite script(s):		14.1
# Prerequisite file(s):		 	N/A
# Description: 					Move the archaea dowloads to a folder containing the archaea genomes
# Output file(s):				


import os, shutil



##########################
# VARIABLES
##########################

# Get the current directory
current_directory = os.getcwd()




# Define the folder that contains the genome downloads
downloads_path = current_directory + "/archaea_genomes/"

# Define the folder to hold the genomes
genome_folder_directory = downloads_path + "archaea_raw_embl/"


print downloads_path
print genome_folder_directory

##########################
# FUNCTIONS
##########################


# Setup new directories
def setup_directory(directory):

	if not os.path.exists(directory):
		os.makedirs(directory)
		print "\nCreated the new directory: %s" % directory
	else:
		print "\nThe directory %s already exists." % directory



# Get a list of the genomes
def get_genomes(directory):
	
	files = []
	
	# For each file in the downloaded genomes
	for eachfile in os.listdir(directory):

		# If the file ends with .embl format
		if eachfile.endswith(".embl"):

			# Append the file to the files list
			files.append(eachfile)

	return files


# Move file to new directory
def move_file(filepath, to_directory, moved_path):

	
	# If the genome file exists
	if os.path.exists(filepath):

		# If the genome file doesnt already exist
		if not os.path.exists(moved_path):
			shutil.move(filepath, to_directory)
		else:
			os.remove(filepath)


def main():

	# Setup the new folder to contain the genome embl files
	setup_directory(genome_folder_directory)

	# Get a list of all the genome files
	genome_files = get_genomes(current_directory)

	#Move the files to the new directory
	print "\nMoving genome files"
	for genome in genome_files:

		genome_path = current_directory + "/" + genome
		moved_path = genome_folder_directory + genome
		move_file(genome_path, genome_folder_directory, moved_path)





# Initiate the script
if __name__ == "__main__":
	main()
