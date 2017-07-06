#!/usr/bin/python

# Script number: 			17.1
# File: 					1 of 1
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2, 4.1, 8.1, 8.2, 8.3, 10.1, 10.2, 10.3, 13.1, 13.2
# Prerequisite file(s):		_filteredGenes.csv, _leaderGenes.csv, _site_4_ratios.csv, _site_4_ratios_t4.csv, SD outputs, CAI outputs
# Description: 				Format data for multivar analysis
# Output file(s):			_all.csv, _reduced.csv, _gene_level.csv


import sys, os, re, shutil, Bio, csv, numpy
import rpy2.robjects as robjects
import multiprocessing as mp
from numpy import array

##########################
# VARIABLES
##########################





genomes_dir = "genome_extractions/cds/"
expression_dir = "outputs/expression/"
expression_genomes_dir = expression_dir + "genomes/"
output_directory = "outputs/multivar/"

codon_map = {"TTT":"F", "TTC":"F", "TTA":"Lt", "TTG":"Lt",
	   "TCT":"St", "TCC":"St", "TCA":"St", "TCG":"St",
	   "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
	   "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
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

aminos = []
for codon in codon_map:
	aminos.append(codon_map[codon])
aminos = sorted(numpy.unique(aminos))



##########################
# FUNCTIONS
##########################

def setupDirectories(direc):

	# If the directory doesnt exist, make new directory
	if not os.path.exists(direc):
		os.makedirs(direc)

# Create a list of all the good genomes
def getGenomes(direc):

	files = []

	# For each file in the downloaded genomes
	for eachfile in os.listdir(direc):

		if eachfile != '.DS_Store':
			# if eachfile in testFile:
				files.append(eachfile[:-4])

	return files, len(files)

# Split the text file containing the CDS for the genome
def splitFile(file, split):
	with open(file) as myfile:
		line = myfile.read()
		genes = line.split(split)

		return genes



def proportion(number, total):

	if total != 0:
		pc = number / float(total)
	else:
		pc = 0
	return pc


def getCAIs(acc):

	CAIvalues = {}


	# # Split the file containing the CAI values and
	# # remove the first and last entry (blank)
	CAIfilePath = expression_genomes_dir + acc + "/" + acc + ".out"
	CAIs = splitFile(CAIfilePath, "\n")
	CAIs.pop()

	CAIline = 0

	for cds in CAIs:

		CAIline +=1

		# For each of the CAIs
		if CAIline > 1:

			# Split the line
			cdsLine = cds.split("\t")

			# Retrieve the gene name and the CAI value
			geneName = cdsLine[0]
			CAIval = cdsLine[8]

			if CAIval != '*****':
				CAIval = float(CAIval)
				CAIvalues[geneName] = CAIval

	return CAIvalues



def get_genome(fourthPath):

	GCs = {}


	with open(fourthPath, "U") as myfile:

		lineNum = 0
		for line in myfile:

			lineNum +=1

			if lineNum >1:
				splits = line.split(",")
				acc = splits[1]
				GC3 = splits[5]

				GCs[acc] = GC3


	return GCs







def get_geneSeq(acc):


	geneSeq = {}
	locus_gene = {}

	genomeFile = genomes_dir + acc + ".txt"

	with open(genomeFile, "U") as myfile:

		line = myfile.read()
		genes = line.split("\n\n")
		genes.pop()

		unknownCount = 0

		# For each of the genes in the list of genes in the genome
		for singleGene in genes:

			# Extract the gene
			gene = re.findall('gene=(.*?);', singleGene)[0]

			# If the gene isnt blank
			if gene != '':

				# Get the gene
				gene = str(gene).replace("['",'').replace("']",'')

			else:

				# Otherwise increase the unknown count
				unknownCount += 1

				# Get the gene known to unknown
				gene = "unknown%s" % unknownCount

			# Get the sequCAIe for the gene
			seq = re.findall('(?<=\n).*', singleGene)[0]

			locus = re.findall('locus_tag=(.*?);', singleGene)[0]

			# print gene, len(seq)

			geneSeq[gene] = seq
			locus_gene[gene] = locus

			# print acc, gene, geneSeq[acc][gene]

	return geneSeq, locus_gene







def get_HE_genes(CAIvals, acc):


	HEgenes = []



	# Get the genes used as highly expressed genes for each genome
	hePath = expression_genomes_dir + acc + "/" + acc + "_he.txt"

	with open(hePath, "U") as myfile:

		line = myfile.read()
		genes = line.split("\n")


		for singleHEgene in genes:

			gene = re.findall('^>(.+)(?=\\t\\t\\t)', singleHEgene)



			if len(gene) != 0:
				HEgenes.append(gene[0])

	return HEgenes



def extract_dg_values(acc):


	dgFileDir = 'outputs/sd_sequences/antiSD_upstream_seq/' + acc +  "/dg_values/"

	genome_sds = {}

	# If there a dG values
	if len(os.listdir(dgFileDir)) != 0:

		dGvalues = []
		fileNo = 0

		# print ("...Extracting dG values for %s" % acc)
		# gene_sequences = get_gene_sequences_for_dg(acc, filteredGenes)


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

					max_dG = [0,0]

					for pos in genome_dGvalues:

						if genome_dGvalues[pos] < max_dG[1]:
							max_dG = [pos,genome_dGvalues[pos]]

					# if max_dG[0] < 30 and max_dG[1] < -3.4535:
					# 	genome_sds[locus] = 1
					# else:
					# 	genome_sds[locus] = 0
					genome_sds[locus] = max_dG[1]

		return genome_sds










def check_exists(acc):

	sd_path = 'outputs/sd_sequences/antiSD_upstream_seq/' + acc + '/'
	cai_path = expression_genomes_dir + acc + "/"


	if os.path.exists(sd_path) and os.path.exists(cai_path):
		return True
	else:
		return False


def get_leaders():

	leaders = {}

	leader_path = 'outputs/leader_genes/_leaderGenes.csv'

	# Open the file containing the filtered genes
	with open(leader_path, "U") as myfile:
		rowNum  = 0

		# For each row in the file
		for row in myfile:

			rowNum += 1

			# Exclude the header row
			if rowNum > 1:

				# Split the row
				leadergene = row.split(',')

				# Get the accession number of the genome
				acc = leadergene[0]

				locus = leadergene[1]


				if acc in leaders:
					leaders[acc].append(locus)
				else:
					leaders[acc] = []
					leaders[acc].append(locus)

	return leaders


def get_gc_content():

	GC = {}
	file_path = 'outputs/ratio_testing/site_4/_site_4_ratios.csv'

	with open(file_path, 'rU') as ratio_file:
		lines = ratio_file.readlines()

		line_count = 0
		for line in lines:
			line_count += 1
			if line_count > 1:
				splits = line.split(',')
				acc = splits[1]
				gc = splits[5]
				GC[acc] = gc

	file_path = 'outputs/ratio_testing/site_4/_site_4_ratios_t4.csv'
	with open(file_path, 'rU') as ratio_file:
		lines = ratio_file.readlines()

		line_count = 0
		for line in lines:
			line_count += 1
			if line_count > 1:
				splits = line.split(',')
				acc = splits[1]
				gc = splits[5]
				GC[acc] = gc

	return GC

def all_genomes(acc, CAIvalues, gene_seqs, genome_sds, HE_gene, locus_gene, genome_leaders):

	genome_stats = {}

	if genome_sds != None:
		for gene in locus_gene:
			if gene in CAIvalues:
				if locus_gene[gene] in genome_sds and gene not in HE_gene:

					sites = [6,7,9,10,12]
					a_content = 0
					t_content = 0

					a_count = 0
					codon_count = 0

					for site in sites:
						if  gene_seqs[gene][site-1] == "A":
							a_content += 1
					a_content = a_content / float(len(sites))


					fourth =0
					# fourth = gene_seqs[gene][3]
					if gene_seqs[gene][3] == "A":
						fourth = 1

					if locus_gene[gene] in genome_leaders:
						leader = 1
					else:
						leader = 0

					for i in range(6, len(gene_seqs[gene])-3, 3):
						codon_count += 1
						if gene_seqs[gene][i] == "A":
							a_count += 1




					genome_stats[locus_gene[gene]] = [fourth, CAIvalues[gene], genome_sds[locus_gene[gene]], a_content, leader, a_count, codon_count]
					# output_line += '%s,%s,%s,%s,%s,%s,%s\n' % (locus_gene[gene], fourth, genome_sds[locus_gene[gene]], CAIvalues[gene], a_content, leader, to_stop)






	return genome_stats


def no_cai_sd(acc, gene_seqs, locus_gene, genome_leaders):

	genomes_no_cai_sd_stats = {}


	for gene in locus_gene:

		sites = [6,7,9,10,12]
		a_content = 0
		t_content = 0

		a_count = 0
		codon_count = 0

		for site in sites:
			if  gene_seqs[gene][site-1] == "A":
				a_content += 1
		a_content = a_content / float(len(sites))


		fourth =0
		# fourth = gene_seqs[gene][3]
		if gene_seqs[gene][3] == "A":
			fourth = 1

		if locus_gene[gene] in genome_leaders:
			leader = 1
		else:
			leader = 0

		for i in range(6, len(gene_seqs[gene])-3, 3):
			codon_count += 1
			if gene_seqs[gene][i] == "A":
				a_count += 1




		genomes_no_cai_sd_stats[locus_gene[gene]] = [fourth, a_content, leader, a_count, codon_count, gene_seqs[gene][3]]


	return genomes_no_cai_sd_stats


def run_genome(acc, fileNumber, leaders,exists):

	print('%s: %s' % (acc, fileNumber))

	CAIvals = {}
	geneSeq = {}

	gene_seqs, locus_gene = get_geneSeq(acc)
	GC = get_gc_content()


	# Get the data if SD exists for the genome
	if exists:
		CAIvalues = getCAIs(acc)
		HE_genes = get_HE_genes(CAIvalues, acc)
		genome_sds = extract_dg_values(acc)
		genome_stats = all_genomes(acc, CAIvalues, gene_seqs, genome_sds, HE_genes, locus_gene, leaders[acc])
	else:
		genome_stats = False

	genomes_no_cai_sd_stats = no_cai_sd(acc, gene_seqs, locus_gene, leaders[acc])




	return acc, genome_stats, genomes_no_cai_sd_stats, exists
	# print_genome(acc, output_line)

def get_table_4_genomes():

	table4_genomes = []

	with open('outputs/gene_filtering/table4genomes.txt', 'rU') as file:

		lines = file.readlines()
		for line in lines:
			line = line.split(',')
			table4_genomes.append(line[1].strip('\n'))

	return(table4_genomes)

def main():


	fileNumber = 0

	setupDirectories(output_directory)
	# output_file = open(output_directory + 'multivar.csv', 'w')
	# #output_file.write('acc,gc,pval,padj,mean_CAI_a,mean_CAI_not_a\n')

	# Get a list of the genomes containing CAI information
	genomes, numGenomes = getGenomes(genomes_dir)

	table4_genomes = get_table_4_genomes()

	leaders = get_leaders()

	numProcessors = (mp.cpu_count()) - 2
	pool = mp.Pool(numProcessors)
	tasks = []

	GC = get_gc_content()

	for acc in genomes:
:
			fileNumber += 1

			exists = check_exists(acc)

			tasks.append((acc,fileNumber,leaders,exists))


	results = []
	for task in tasks:
		results.append(pool.apply_async(run_genome, task))


	output_file_2 = open(output_directory + '_all.csv', 'w')
	output_file_2.write('acc,fourth,cai,sd,a_content,leader,gc,a_codons,trans_table\n')

	output_file_3 = open(output_directory + '_reduced.csv', 'w')
	output_file_3.write('acc,fourth,a_content,leader,gc,a_codons,trans_table,a_count,non_a_count\n')


	output_file_4 = open(output_directory + '_gene_level.csv', 'w')
	output_file_4.write('acc,locus,fourth,a_content,local_a_content,leader,gc,a_codons,trans_table\n')

	output_file_5 = open(output_directory + '_gene_level2.csv', 'w')
	output_file_5.write('acc,locus,fourth,a_content,local_a_content,leader,gc,a_codons,trans_table\n')

	for result in results:
		acc,  genome_stats, genomes_no_cai_sd_stats, exists = result.get()
		# print_genome(acc, output_line)

		if exists:
			if len(genome_stats) > 0:

				fourth = []
				cai = []
				sd = []
				a_content = []
				leaders = []

				a_total = 0
				codon_total = 0

				for locus in genome_stats:

					a_total += genome_stats[locus][5]
					codon_total += genome_stats[locus][6]

					fourth.append(genome_stats[locus][0])
					cai.append(genome_stats[locus][1])

					if genome_stats[locus][2] < -3.4535:
						sd.append(1)
					else:
						sd.append(0)

					a_content.append(genome_stats[locus][3])
					leaders.append(genome_stats[locus][4])




				m_fourth = numpy.mean(fourth)
				m_cai = numpy.mean(cai)
				m_sd = numpy.mean(sd)
				m_a_content = numpy.mean(a_content)
				m_leaders = numpy.mean(leaders)

				a_codons = a_total / float(codon_total)

				if acc in table4_genomes:
					table = 'table4'
				else:
					table = 'table11'


				output_line = '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (acc, m_fourth, m_cai, m_sd, m_a_content, m_leaders, GC[acc], a_codons, table)
				output_file_2.write(output_line)





		fourth = []
		cai = []
		sd = []
		a_content = []
		leaders = []

		a_total = 0
		codon_total = 0

		fourth_non_a_count = 0
		fourth_a_count = 0

		for locus in genomes_no_cai_sd_stats:


			a_total += genomes_no_cai_sd_stats[locus][3]

			if genomes_no_cai_sd_stats[locus][0] != 1:
				fourth_non_a_count += 1
			else:
				fourth_a_count += 1

			codon_total += genomes_no_cai_sd_stats[locus][4]

			fourth.append(genomes_no_cai_sd_stats[locus][0])


			a_content.append(genomes_no_cai_sd_stats[locus][1])
			leaders.append(genomes_no_cai_sd_stats[locus][2])

		m_fourth = numpy.mean(fourth)

		m_a_content = numpy.mean(a_content)
		m_leaders = numpy.mean(leaders)

		a_codons = a_total / float(codon_total)

		if acc in table4_genomes:
			table = 'table4'
		else:
			table = 'table11'




		output_line = '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (acc, m_fourth, m_a_content, m_leaders, GC[acc], a_codons, table, fourth_a_count, fourth_non_a_count)
		output_file_3.write(output_line)


		for locus in genomes_no_cai_sd_stats:

			fourth_site = genomes_no_cai_sd_stats[locus][0]
			a_content = genomes_no_cai_sd_stats[locus][1]
			leader = genomes_no_cai_sd_stats[locus][2]



			output_line = '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (acc, locus, fourth_site, m_a_content, a_content, leader, GC[acc], a_codons, table)
			output_file_4.write(output_line)

			output_line = '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (acc, locus, genomes_no_cai_sd_stats[locus][5], m_a_content, a_content, leader, GC[acc], a_codons, table)
			output_file_5.write(output_line)


	output_file_2.close()
	output_file_3.close()
	output_file_4.close()




#################################

# Initiate the filtering
if __name__ == "__main__":
	main()
