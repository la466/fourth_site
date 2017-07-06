#!/usr/bin/python

# Script number: 			12.1
# File: 					1 of 1
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv, table4genomes.txt
# Description: 				Get codon frequencies and see how common they are after methionine
# Output file(s):			_met_codon_pair.csv, _met_codon_pair_summary.csv

import csv, re, numpy, os
import rpy2.robjects as robjects

############
# VARIABLES
############

filtered_genes_file_path = "outputs/gene_filtering/_filteredGenes.csv"
genomes_cds_path = "genome_extractions/cds/"
table4_genomes_path = "outputs/gene_filtering/table4genomes.txt"

output_dir = "outputs/met_pair/"
output_file_path = output_dir + "_met_codon_pair.csv"
output_file2_path = output_dir + "_met_codon_pair_summary.csv"

codons = [
		"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
		"CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
		"GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
		"TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT",
	]

codon_map = {
		"TTT":"F", "TTC":"F", "TTA":"Lt", "TTG":"Lt",
      	"TCT":"St", "TCC":"St", "TCA":"St", "TCG":"St",
       	"TAT":"Y", "TAC":"Y", "TAA":"Stop", "TAG":"Stop",
       	"TGT":"C", "TGC":"C", "TGA":"Stop", "TGG":"W",
       	"CTT":"Lc", "CTC":"Lc", "CTA":"Lc", "CTG":"Lc",
       	"CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       	"CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       	"CGT":"Rc", "CGC":"Rc", "CGA":"Rc", "CGG":"Rc",
       	"ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       	"ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       	"AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       	"AGT":"Sa", "AGC":"Sa", "AGA":"Ra", "AGG":"Ra",
       	"GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       	"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       	"GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       	"GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
 	}

aminos = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "Lc", "Lt", "M", "N", "P", "Q", "Ra", "Rc", "Sa", "St", "Stop", "T", "V", "W", "Y"]


testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750"]


###########
# FUNCTIONS
###########

# Setup new directories
def setupDirectory(directory):

	if not os.path.exists(directory):
		os.makedirs(directory)
		print "Created the new directory: %s" % directory

def get_filtered_genes(file_path):

	filteredGenes = {}

	with open(file_path, 'U') as myfile:
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


def get_table4_genomes(file_path):

	with open(file_path, "U") as myfile:

		table4_genomes = []

		lines = myfile.readlines()
		for line in lines:
			line = line.replace("\n", "")
			table4_acc = line.split(",")[1]

			table4_genomes.append(table4_acc)

		return table4_genomes



def read_genome(acc, filtered_genes, table4_genomes, genomes_cds_path):

	cds = []


	genome_cds_path = genomes_cds_path + acc + ".txt"

	with open(genome_cds_path, "U") as genome_file:



		genome_cds = genome_file.read()
		all_cds = genome_cds.split('\n\n')
		all_cds.pop()

		for single_cds in all_cds:

			# Get the locus of the cds
			cds_locus = re.findall('locus_tag=(.+?);', single_cds)[0]

			# Make sure the cds has passed filtering
			if cds_locus in filtered_genes[acc]:

				seq = re.findall('\n(.+)', single_cds)[0]

				cds.append(seq)

	return cds


def codon_repeat_analysis(acc, cds, codons, output_file, chi_vals, chi_vals_amino):

	codon_counts = {}
	codon_repeats = {}

	for codon in codons:
		codon_counts[codon] = 0
		codon_repeats[codon] = 0

	codon_count = 0
	codon_repeat_count = 0


	amino_repeats = {}
	amino_counts = {}
	for amino in aminos:
		amino_counts[amino] = 0
		amino_repeats[amino] = 0

	amino_count = 0
	amino_repeat_count = 0

	for single_cds in cds:


		# Get query sequence (not including start codon, second codon or stop)
		query_seq = single_cds[6:-3]
		codon_count += len(query_seq)/3
		amino_count += len(query_seq)/3

		# Get each codon in the query sequence
		for i in range(0,len(query_seq),3):
			codon = query_seq[i:i+3]
			codon_counts[codon] += 1

			amino_counts[codon_map[codon]] += 1


			# If the codon is not the last codon
			if (i+3) < len(query_seq):

				# Add to the number of pairs
				codon_repeat_count += 1
				amino_repeat_count += 1

				# If the current codon is methionine,
				# look for the next codon
				if codon == "ATG":
					next_codon = query_seq[i+3:i+6]
					codon_repeats[next_codon] += 1
					amino_repeats[codon_map[next_codon]] += 1





	file_string = "%s" % acc


	met_prop = codon_counts["ATG"] / float(codon_count)




	for codon in sorted(codon_counts.iterkeys()):

		# Genome codon proportion
		if codon_counts[codon] != 0:
			codon_prop = codon_counts[codon] / float(codon_count)
		else:
			codon_prop = 0

		# Genome expected met-codon pair proportion
		codon_met_exp_prop = met_prop * codon_prop

		# Get expected number of met-codon pairs
		codon_met_exp_count = codon_met_exp_prop * float(codon_repeat_count)

		# print codon, codon_prop, codon_met_exp_count

		if codon_met_exp_count != 0:
			chi_square_val = ((codon_repeats[codon] - codon_met_exp_count)**2) / float(codon_met_exp_count)
		else:
			chi_square_val = 0

		codon_diff = codon_repeats[codon] - codon_met_exp_count
		chi_vals[codon][0].append(codon_diff)
		chi_vals[codon][1].append(chi_square_val)

		file_string += ",%s,%s,%s" % (codon_met_exp_count, codon_repeats[codon], chi_square_val)

	file_string += "\n"
	output_file.write(file_string)


	for amino in sorted(amino_counts.iterkeys()):
		if amino_counts[amino] != 0:
			amino_prop = amino_counts[amino] / float(amino_count)
		else:
			amino_prop = 0


		# Genome expected met-amino pair proportion
		amino_met_exp_prop = met_prop * amino_prop


		# Get expected number of met-amino pairs
		amino_met_exp_count = amino_met_exp_prop * float(amino_repeat_count)


		# print amino, amino_prop, amino_met_exp_count

		if amino_met_exp_count != 0:
			chi_square_val = ((amino_repeats[amino] - amino_met_exp_count)**2) / float(amino_met_exp_count)
		else:
			chi_square_val = 0


		amino_diff = amino_repeats[amino] - amino_met_exp_count

		chi_vals_amino[amino][0].append(amino_diff)
		chi_vals_amino[amino][1].append(chi_square_val)



	return chi_vals, chi_vals_amino


def calc_chi_vals(chi_vals, chi_vals_amino, output_file, output_file2):

	codon_collate = {}
	for codon in codons:
		codon_collate[codon] = []

	file_string = "mean_obs-exp"
	for codon in sorted(chi_vals.iterkeys()):
		file_string += ",,,%s" % numpy.mean(chi_vals[codon][0])
		codon_collate[codon].append(numpy.mean(chi_vals[codon][0]))
	file_string += "\n"
	output_file.write(file_string)


	file_string = "chisq_sum"
	for codon in sorted(chi_vals.iterkeys()):
		file_string += ",,,%s" % sum(chi_vals[codon][1])
		codon_collate[codon].append(sum(chi_vals[codon][1]))
	file_string += "\n"
	output_file.write(file_string)

	file_string = "p_val"
	for codon in sorted(chi_vals.iterkeys()):
		chi_sum = sum(chi_vals[codon][1])
		p_value = robjects.r('pchisq(%s,%s,lower.tail=FALSE)' % (chi_sum, len(chi_vals[codon][1])-1))[0]
		file_string += ",,,%s" % p_value
		codon_collate[codon].append(p_value)
	output_file.write(file_string)


	for codon in sorted(codon_collate.iterkeys()):
		file_string = "%s,%s,%s,%s,,%s\n" % (codon,codon_collate[codon][0],codon_collate[codon][1],codon_collate[codon][2],codon_map[codon])
		output_file2.write(file_string)
	output_file2.write("\n")


	amino_collate = {}
	for amino in aminos:
		amino_collate[amino] = []

	for amino in sorted(chi_vals_amino.iterkeys()):

		amino_collate[amino].append(numpy.mean(chi_vals_amino[amino][0]))

	for amino in sorted(chi_vals_amino.iterkeys()):
		amino_collate[amino].append(sum(chi_vals_amino[amino][1]))

	for amino in sorted(chi_vals_amino.iterkeys()):
		chi_sum = sum(chi_vals_amino[amino][1])
		p_value = robjects.r('pchisq(%s,%s,lower.tail=FALSE)' % (chi_sum, len(chi_vals_amino[amino][1])-1))[0]
		amino_collate[amino].append(p_value)


	for amino in sorted(amino_collate.iterkeys()):
		file_string = "%s,%s,%s,%s\n" % (amino,amino_collate[amino][0],amino_collate[amino][1],amino_collate[amino][2])
		output_file2.write(file_string)



def main():


	setupDirectory(output_dir)

	# Get a list of the filtered genes
	filtered_genes = get_filtered_genes(filtered_genes_file_path)

	# Get a list of table 4 genomes
	table4_genomes = get_table4_genomes(table4_genomes_path)

	output_file = open(output_file_path, "w")
	output_file.write("acc")
	for codon in codons:
		file_string = ",%s_exp,%s_obs,%s_chisq_calc(O-E^2/E)" % (codon,codon,codon)
		output_file.write(file_string)
	output_file.write("\n")

	output_file2 = open(output_file2_path, "w")
	output_file2.write("codon,mean_obs-exp,chisq_sum,p_value,,amino\n")

	genome_count = 0

	chi_vals = {}
	for codon in codons:
		chi_vals[codon] = [[],[]]

	chi_vals_amino = {}
	for amino in aminos:
		chi_vals_amino[amino] = [[],[]]

	for acc in filtered_genes:

		if acc not in table4_genomes:
		# if acc == "AP009152":
		# if acc in testFiles:
			genome_count+=1
			print "-" * 20
			print "Genome %d of %d\n%s\n" % (genome_count, len(filtered_genes), acc)

			# Extract the cds
			cds = read_genome(acc, filtered_genes, table4_genomes, genomes_cds_path)

			# Get codon frequency
			chi_vals, chi_vals_amino = codon_repeat_analysis(acc, cds, codons, output_file, chi_vals, chi_vals_amino)

	calc_chi_vals(chi_vals, chi_vals_amino, output_file, output_file2)


if __name__ == "__main__":

	main()
