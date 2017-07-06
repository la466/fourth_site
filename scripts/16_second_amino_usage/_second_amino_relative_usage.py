#!/usr/bin/python

# Script number: 			16.1
# File: 					1 of 2
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv, table4genomes.txt
# Description: 				Look at the usage of amino acids in the second position
# Output file(s):			_second_amino_usage.csv, _second_amino_usage_summary.csv


import csv, re, numpy, os
import rpy2.robjects as robjects

############
# Variables
############

filtered_genes_file_path = "outputs/gene_filtering/_filteredGenes.csv"
genomes_cds_path = "genome_extractions/cds/"
table4_genomes_path = "outputs/gene_filtering/table4genomes.txt"

output_dir = "outputs/second_amino/"
output_file_path = output_dir + "_second_amino_usage.csv"
output_file2_path = output_dir + "second_amino_usage_summary.csv"


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
# Functions
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


def codon_repeat_analysis(acc, cds, codons, output_file, chi_vals):

	amino_counts = {}
	second_amino_counts = {}

	for codon in codons:
		amino_counts[codon_map[codon]] = 0
		second_amino_counts[codon_map[codon]] = 0

	amino_count = 0
	cds_count = 0

	# Get the number of each amino acid used
	for single_cds in cds:

		cds_count += 1

		# Get query sequence (not including start codon, second codon or stop)
		query_seq = single_cds[3:-3]
		amino_count += len(query_seq)/3

		second_amino_counts[codon_map[query_seq[0:3]]] += 1

		# Get each codon in the query sequence
		for i in range(0,len(query_seq),3):
			codon = query_seq[i:i+3]
			amino_counts[codon_map[codon]] += 1





	file_string = "%s" % acc






	for amino in sorted(amino_counts.iterkeys()):

		# Genome amino proportion
		if amino_counts[amino] != 0:
			amino_prop = amino_counts[amino] / float(amino_count)
		else:
			amino_prop = 0


		amino_exp = amino_prop * cds_count
		amino_obs = second_amino_counts[amino]

		if amino_exp != 0:
			chi_square_val = ((amino_obs - amino_exp)**2 / amino_exp)
		else:
			chi_square_val = 0


		chi_vals[amino][0].append(amino_obs - amino_exp)
		chi_vals[amino][1].append(chi_square_val)



		file_string += ",%s,%s,%s" % (amino_obs, amino_exp, chi_square_val)

	serine_ratio = amino_counts["Sa"]/float(amino_counts["St"])
	serine_second_ratio = second_amino_counts["Sa"]/float(second_amino_counts["St"])

	# proportion of A starting serine out of total serine
	sera_ratio = (amino_counts["Sa"]/float(amino_counts["Sa"] + amino_counts["St"]))
	sert_ratio = (amino_counts["St"]/float(amino_counts["Sa"] + amino_counts["St"]))

	# Expected A starting serine: number of serine at second site * serine ratio
	exp_serinea_second = ((second_amino_counts["Sa"] + second_amino_counts["St"]) * sera_ratio)
	exp_serinet_second = ((second_amino_counts["Sa"] + second_amino_counts["St"]) * sert_ratio)
	# print exp_serinea_second

	# Observed serine at the second position
	obs_serinea_second = second_amino_counts["Sa"]
	obs_serinet_second = second_amino_counts["St"]
	# print obs_serinea_second

	obs_exp_sera = ((obs_serinea_second - exp_serinea_second)**2) / float(exp_serinea_second)
	obs_exp_sert = ((obs_serinet_second - exp_serinet_second)**2) / float(exp_serinet_second)

	chi_vals["sera"][0].append(obs_serinea_second - exp_serinea_second)
	chi_vals["sera"][1].append(obs_exp_sera)

	chi_vals["sert"][0].append(obs_serinet_second - exp_serinet_second)
	chi_vals["sert"][1].append(obs_exp_sert)

	arginine_ratio = amino_counts["Ra"]/float(amino_counts["Rc"])
	arginine_second_ratio = second_amino_counts["Ra"]/float(second_amino_counts["Rc"])

	# proportion of A starting arginine out of total arginine
	arga_ratio = (amino_counts["Ra"]/float(amino_counts["Ra"] + amino_counts["Rc"]))
	argc_ratio = (amino_counts["Rc"]/float(amino_counts["Ra"] + amino_counts["Rc"]))

	# Expected A starting arginine: number of arginine at second site * arginine ratio
	exp_argininea_second = ((second_amino_counts["Sa"] + second_amino_counts["Rc"]) * arga_ratio)
	exp_argininec_second = ((second_amino_counts["Sa"] + second_amino_counts["Rc"]) * argc_ratio)
	# print exp_argininea_second

	# Obargved arginine at the second position
	obs_argininea_second = second_amino_counts["Ra"]
	obs_argininec_second = second_amino_counts["Rc"]
	# print obs_argininea_second

	obs_exp_arga = ((obs_argininea_second - exp_argininea_second)**2) / float(exp_argininea_second)
	obs_exp_argc = ((obs_argininec_second - exp_argininec_second)**2) / float(exp_argininec_second)

	chi_vals["arga"][0].append(obs_argininea_second - exp_argininea_second)
	chi_vals["arga"][1].append(obs_exp_arga)

	chi_vals["argc"][0].append(obs_argininec_second - exp_argininec_second)
	chi_vals["argc"][1].append(obs_exp_argc)

	file_string += ",%s,%s,%s,%s,%s,%s,%s,%s" % (serine_ratio, serine_second_ratio,obs_serinea_second, exp_serinea_second,obs_exp_sera,obs_serinet_second, exp_serinet_second,obs_exp_sert)
	file_string += ",%s,%s,%s,%s,%s,%s,%s,%s\n" % (arginine_ratio, arginine_second_ratio,obs_argininea_second, exp_argininea_second,obs_exp_arga,obs_argininec_second, exp_argininec_second,obs_exp_sert)

	output_file.write(file_string)

	# print "Total serine: %s" % (amino_counts["Sa"] + amino_counts["St"])
	# print "Total A serine: %s" % amino_counts["Sa"]
	# print "A serine prop: %s\n" % sera_ratio

	# print "Total second serine: %s" % (second_amino_counts["Sa"] + second_amino_counts["St"])
	# print "Obs second A serine: %s" % second_amino_counts["Sa"]
	# print "Expected Second A serine: %s" % ((second_amino_counts["Sa"] + second_amino_counts["St"]) * sera_ratio)

	return chi_vals


def calc_chi_vals(chi_vals, output_file, output_file2):

	sera_mean_obs_exp = numpy.mean(chi_vals["sera"][0])
	sera_chi_sum = sum(chi_vals["sera"][1])
	sera_p_value = robjects.r('pchisq(%s,%s,lower.tail=FALSE)' % (sera_chi_sum, len(chi_vals["sera"][1])-1))[0]

	sert_mean_obs_exp = numpy.mean(chi_vals["sert"][0])
	sert_chi_sum = sum(chi_vals["sert"][1])
	sert_p_value = robjects.r('pchisq(%s,%s,lower.tail=FALSE)' % (sert_chi_sum, len(chi_vals["sert"][1])-1))[0]


	arga_mean_obs_exp = numpy.mean(chi_vals["arga"][0])
	arga_chi_sum = sum(chi_vals["arga"][1])
	arga_p_value = robjects.r('pchisq(%s,%s,lower.tail=FALSE)' % (arga_chi_sum, len(chi_vals["arga"][1])-1))[0]

	argc_mean_obs_exp = numpy.mean(chi_vals["argc"][0])
	argc_chi_sum = sum(chi_vals["argc"][1])
	argc_p_value = robjects.r('pchisq(%s,%s,lower.tail=FALSE)' % (argc_chi_sum, len(chi_vals["argc"][1])-1))[0]


	amino_block_aminos = ["sera", "sert", "arga", "argc"]

	amino_collate = {}
	for codon in codons:
		amino_collate[codon_map[codon]] = []

	file_string = "mean_obs-exp"
	for amino in sorted(chi_vals.iterkeys()):
		if amino not in amino_block_aminos:
			file_string += ",,,%s" % numpy.mean(chi_vals[amino][0])
			amino_collate[amino].append(numpy.mean(chi_vals[amino][0]))
	file_string += ",,,,,%s,,,%s" % (sera_mean_obs_exp, sert_mean_obs_exp)
	file_string += ",,,,,%s,,,%s" % (arga_mean_obs_exp, argc_mean_obs_exp)
	file_string += "\n"
	output_file.write(file_string)


	file_string = "chisq_sum"
	for amino in sorted(chi_vals.iterkeys()):
		if amino not in amino_block_aminos:
			file_string += ",,,%s" % sum(chi_vals[amino][1])
			amino_collate[amino].append(sum(chi_vals[amino][1]))
	file_string += ",,,,,%s,,,%s" % (sera_chi_sum, sert_chi_sum)
	file_string += ",,,,,%s,,,%s" % (arga_chi_sum, argc_chi_sum)
	file_string += "\n"
	output_file.write(file_string)

	file_string = "p_val"
	for amino in sorted(chi_vals.iterkeys()):
		if amino not in amino_block_aminos:
			chi_sum = sum(chi_vals[amino][1])
			p_value = robjects.r('pchisq(%s,%s,lower.tail=FALSE)' % (chi_sum, len(chi_vals[amino][1])-1))[0]
			file_string += ",,,%s" % p_value
			amino_collate[amino].append(p_value)
	file_string += ",,,,,%s,,,%s" % (sera_p_value, sert_p_value)
	file_string += ",,,,,%s,,,%s" % (arga_p_value, argc_p_value)
	output_file.write(file_string)


	for amino in sorted(amino_collate.iterkeys()):
		if amino not in amino_block_aminos:
			file_string = "%s,%s,%s,%s\n" % (amino,amino_collate[amino][0],amino_collate[amino][1],amino_collate[amino][2])
			output_file2.write(file_string)
	ser_string = "\nser_a_from_genome_serine_prop,%s,%s,%s\n" % (sera_mean_obs_exp, sera_chi_sum, sera_p_value)
	ser_string += "ser_t_from_genome_serine_prop,%s,%s,%s\n" % (sert_mean_obs_exp, sert_chi_sum, sert_p_value)
	arg_string = "\narg_a_from_genome_arg_prop,%s,%s,%s\n" % (arga_mean_obs_exp, arga_chi_sum, arga_p_value)
	arg_string += "arg_c_from_genome_arg_prop,%s,%s,%s\n" % (argc_mean_obs_exp, argc_chi_sum, argc_p_value)
	output_file2.write(ser_string)
	output_file2.write(arg_string)

def main():

	setupDirectory(output_dir)

	# Get a list of the filtered genes
	filtered_genes = get_filtered_genes(filtered_genes_file_path)

	# Get a list of table 4 genomes
	table4_genomes = get_table4_genomes(table4_genomes_path)

	output_file = open(output_file_path, "w")
	output_file.write("acc")
	for amino in aminos:
		file_string = ",%s_obs,%s_exp,%s_(O-E)^2/E" % (amino,amino,amino)
		output_file.write(file_string)
	output_file.write(",serine_ratio(A/T),serine_second_ratio(A/T),obs_second_sera,exp_second_sera_from_serine_usage,sera_(O-E)^2/E,obs_second_sert,exp_second_sert_from_serine_usage,sert_(O-E)^2/E")
	output_file.write(",arginine_ratio(A/C),arginine_second_ratio(A/C),obs_second_arga,exp_second_arga_from_arginine_usage,arga_(O-E)^2/E,obs_second_argc,exp_second_argc_from_arginine_usage,argc_(O-E)^2/E\n")

	output_file2 = open(output_file2_path, "w")
	output_file2.write("amino,mean_obs-exp,chisq_sum,p_value\n")



	genome_count = 0

	chi_vals = {}
	for codon in codons:
		chi_vals[codon_map[codon]] = [[],[]]

	chi_vals["sera"] = [[],[]]
	chi_vals["sert"] = [[],[]]
	chi_vals["arga"] = [[],[]]
	chi_vals["argc"] = [[],[]]

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
			chi_vals = codon_repeat_analysis(acc, cds, codons, output_file, chi_vals)

	calc_chi_vals(chi_vals, output_file, output_file2)


if __name__ == "__main__":

	main()
