#!/usr/bin/python

# Script number: 			9.1
# File: 					1 of 1
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Comapre the changes from each nucleotide in the 5' domain between E coli and Shigella
# Output file(s):			_compare_ecoli_shigella.csv


import os, re, shutil, csv



##########################
# VARIABLES
##########################

ecoli_acc = "AE005174"
shigella_acc = "AE005674"
output_directory = "outputs/species_compare/"
genomes_directory = "genome_extractions/cds/"


##########################
# FUNCTIONS
##########################

# Setup new directories
def setupDirectory(directory):

	if not os.path.exists(directory):
		os.makedirs(directory)
		print "Created the new directory: %s" % directory


# Get a list of the filtered genes
def getFilteredGenes():

	filteredGenes = {}

	filteredGenesFile = "outputs/gene_filtering/_filteredGenes.csv"
	with open(filteredGenesFile, 'U') as myfile:
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
	

# Rewrite the genome file to the correct fasta format for blast
def genome_to_fasta(acc, genome_file, blast_directory, gene_dictionary):
	

	fasta_file_path = blast_directory + acc + "/" + acc + ".fa"
	fasta_file = open(fasta_file_path, "w")


	with open(genome_file, "U") as myfile:

		no_gene = 0
		ids = []

		lines = myfile.read()
		genes = re.findall('>.*\n.*(?=\n\n)', lines)

		for gene in genes:

			gene_info = re.findall('>.*(?=\n)', gene)[0]
			# print gene_info + "\n\n"

			#gene_info_splits = gene_info.split(";")
			#print gene_info_splits

			acc = re.findall('(?<=acc=).*(?=;locus)', gene_info)[0]
			locus = re.findall('(?<=locus_tag=).*(?=;gene=)', gene_info)[0]

			


			gene_name = re.findall("(?<=gene=\[\').*(?=\'\])", gene_info)
			#print locus
			if gene_name:
				gene_name = gene_name[0]
			else:
				no_gene += 1
				gene_name = "no_gene_name_%d" % no_gene
			
			seq = re.findall('(?<=\n).*$', gene)[0]

			fileLine = ">%s|%s\n%s\n" % (acc, locus, seq)
			fasta_file.write(fileLine)

			gene_dictionary[locus] = seq

	fasta_file.close()


# Get the orthologs from the output file
def parse_blast_output(blast_output_file):

	blast_orthologs = {}


	with open(blast_output_file, "U") as myfile:

		lines = myfile.readlines()
		
		for match in lines:
			match = match.strip("\n")

			splits = match.split(",")

			shigella_locus = re.findall("(?<=AE005674\|).*",splits[0])[0]
			ecoli_locus = re.findall("(?<=AE005174\|).*",splits[1])[0]
			match_percent = splits[2]
			e_val = splits[10]
			score = splits[11]

			if shigella_locus in blast_orthologs:
				if e_val >= blast_orthologs[shigella_locus][2]:
					if score > blast_orthologs[shigella_locus][3]:
						blast_orthologs[shigella_locus] = [ecoli_locus, match_percent, e_val, score]
			else:
				blast_orthologs[shigella_locus] = [ecoli_locus, match_percent, e_val, score]

	return blast_orthologs


def compare_orthologs(blast_orthologs, ecoli_genes, shigella_genes, shigella_acc, ecoli_acc, filteredGenes):

	changes = {}

	nonSynchanges = {}
	changeTos = {}

	bases = ["A", "C", "T", "G"]
	positions = range(0,33,3)

	for base in bases:
		changes[base] = []

	for position in positions:
		nonSynchanges[position] = {}
		changeTos[position] = {}
		for base in bases:
			nonSynchanges[position][base] = []
			changeTos[position][base] = {}
			for base2 in bases:
				changeTos[position][base][base2] = []
	
	for locus in blast_orthologs:

		shigella_locus = locus
		ecoli_locus = blast_orthologs[locus][0]

		if shigella_locus in shigella_genes and ecoli_locus in ecoli_genes:
			
			if shigella_locus in filteredGenes[shigella_acc] and ecoli_locus in filteredGenes[ecoli_acc]:

				shigella_seq = shigella_genes[shigella_locus]
				ecoli_seq = ecoli_genes[ecoli_locus]

				if shigella_seq[3] != ecoli_seq[3]:
					if ecoli_seq[3] in bases:
						changes[ecoli_seq[3]].append(1)
				else:
					if ecoli_seq[3] in bases:
						changes[ecoli_seq[3]].append(0)


				for i in positions:
					if ecoli_seq[i] == ecoli_seq[3]:
						if shigella_seq[i] != ecoli_seq[i]:
							if ecoli_seq[3] in bases:
								nonSynchanges[i][ecoli_seq[3]].append(1)
						else:
							if ecoli_seq[3] in bases:
								nonSynchanges[i][ecoli_seq[3]].append(0)

				
				for i in positions:
					if shigella_seq[i] != ecoli_seq[i]:
						changeTos[i][ecoli_seq[i]][shigella_seq[i]].append(1)
					


	for pos in sorted(changeTos):
		for ecoli_base in changeTos[pos]:
			for shigella_base in changeTos[pos][ecoli_base]:
				# print pos, ecoli_base, shigella_base, len(changeTos[pos][ecoli_base][shigella_base])
				pass

	return changes, nonSynchanges, changeTos

#################################

def main():

	# Setup the directory to contain the blast analysis
	blast_directory = output_directory + "blast/"
	setupDirectory(blast_directory)

	# Setup the directories for each genome
	ecoli_directory = blast_directory + ecoli_acc
	shigella_directory = blast_directory + shigella_acc

	setupDirectory(ecoli_directory)
	setupDirectory(shigella_directory)

	# Get the filtered genes
	filtered_genes = getFilteredGenes()

	# Copy the fasta-format genomes to the new directories
	ecoli_genome = genomes_directory + ecoli_acc + ".txt"
	shigella_genome = genomes_directory + shigella_acc + ".txt"

	shutil.copy(ecoli_genome, ecoli_directory)
	shutil.copy(shigella_genome, shigella_directory)

	# Rewrite the genomes to correct fasta formats and get new file paths
	ecoli_genes = {}
	shigella_genes = {}
	genome_to_fasta(ecoli_acc, ecoli_genome, blast_directory, ecoli_genes)
	genome_to_fasta(shigella_acc, shigella_genome, blast_directory, shigella_genes)

	ecoli_blast_fasta = blast_directory + ecoli_acc + "/" + ecoli_acc + ".fa"
	shigella_blast_fasta = blast_directory + shigella_acc + "/" + shigella_acc + ".fa"


	# Create blast database of the ecoli genome
	#shellCommand = "makeblastdb -in %s -parse_seqids -dbtype nucl" % (ecoli_blast_fasta)
	shellCommand = "makeblastdb -in %s -dbtype nucl" % (ecoli_blast_fasta)

	os.system(shellCommand)




	# Blast the shigella genome against the ecoli genome
	print "\n...Blasting the two genomes\n"
	blast_output = blast_directory + "blast_results.out"
	shellCommand = "blastn -db %s -query %s -out %s -outfmt 10" % ("outputs/species_compare/blast/AE005174/AE005174.fa", shigella_blast_fasta, blast_output)
	os.system(shellCommand)

	blast_orthologs = parse_blast_output(blast_output)
	changes, nonSynchanges, changeTos = compare_orthologs(blast_orthologs, ecoli_genes, shigella_genes, shigella_acc, ecoli_acc, filtered_genes)

	
	# for base in changes:
	# 	print base, sum(changes[base]) / float(len(changes[base]))

	
	compareFilePath = output_directory + "_compare_ecoli_shigella.csv"
	compareFile = open(compareFilePath, "w")
	compareFile.write("pos,A,C,T,G\n")

	for pos in sorted(nonSynchanges):

		fileLine = "%s" % (pos+1)
		for base in nonSynchanges[pos]:
			if len(nonSynchanges[pos][base]) == 0:
				changes = 0
			else:
				changes = sum(nonSynchanges[pos][base]) / float(len(nonSynchanges[pos][base]))
			
			fileLine += ",%s" % changes

		fileLine += "\n"

		compareFile.write(fileLine)

	compareFile.close()



#################################

# Initiate
if __name__ == "__main__":
	main()

