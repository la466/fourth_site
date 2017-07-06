#!/usr/bin/python

# Script number: 			11.1
# File: 					1 of 1
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv, table4genomes.txt
# Description: 				Look at whether any presence of a gene upstream is affecting fourth site A content
# Output file(s):			_upstream_cds.csv, _upstream_cds_distances.csv



import csv, re, numpy, os
import rpy2.robjects as robjects


############
# VARIABLES
############

filtered_genes_file_path = "outputs/gene_filtering/_filteredGenes.csv"
genomes_cds_path = "genome_extractions/cds/"
table4_genomes_path = "outputs/gene_filtering/table4genomes.txt"


output_directory = "outputs/upstream_cds/"
output_file_path = output_directory + "_upstream_cds.csv"
output2_file_path = output_directory + "_upstream_cds_distance.csv"

testFiles = ["AE000511", "AE006914", "AM711867", "CP000030","CP000107","CP000232","CP000237","CP000750","CP000925","CP001618","CP001628","CP001643","CP001819","CP001821","CP001867","CP001964","CP002593","CP002810","FO117623","FO203431","HE983995"]

distances = [10,20,30,40,50,100,150,200]

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

	cds = {}


	genome_cds_path = genomes_cds_path + acc + ".txt"

	with open(genome_cds_path, "U") as genome_file:



		genome_cds = genome_file.read()
		all_cds = genome_cds.split('\n\n')
		all_cds.pop()

		cds_num = 0

		for single_cds in all_cds:

			cds_num +=1

			# Get the locus of the cds
			cds_locus = re.findall('locus_tag=(.+?);', single_cds)[0]
			cds_location = re.findall('loc=(.+?);', single_cds)[0]
			cds_start = re.findall('^(\d+)(?=..)', cds_location)[0]
			cds_end = re.findall('(?<=..)\d+$', cds_location)[0]
			cds_strand = re.findall('complement=(.+?);', single_cds)[0]
			cds_seq = re.findall('\n(.+)', single_cds)[0]



			cds[cds_num] = [cds_locus,cds_start,cds_end,cds_strand,cds_seq]



	return cds


def upstream_analysis_all(acc, cds, filtered_genes, output_file, output_file2, dist_a_prop_lead, dist_a_prop_lag, chi_vals, obs_exp_diff):

	lead_strand = {}
	lag_strand = {}
	lead_count = 0
	lag_count = 0

	cds_same = []
	cds_opposite = []

	same_a_dist = []
	same_nota_dist = []
	opp_a_dist =[]
	opp_nota_dist = []


	obs_cds_with_a_with_upstream_cds_same_strand = 0
	total_with_cds_upstream = 0
	obs_cds_with_a_with_upstream_cds = 0
	total_with_cds_upstream_on_same_strand = 0

	for cds_id in cds:

		cds_locus = cds[cds_id][0]
		cds_start = cds[cds_id][1]
		cds_end = cds[cds_id][2]
		cds_strand = cds[cds_id][3]
		cds_seq = cds[cds_id][4]




		if cds_strand == "0":
			lead_count += 1
			lead_strand[lead_count] = [cds_locus,cds_start,cds_end,cds_strand,cds_seq]
		elif cds_strand == "1":
			lag_count += 1
			lag_strand[lag_count] = [cds_locus,cds_start,cds_end,cds_strand,cds_seq]


		# Look at the prop A of non-overlapping genes (either strand)
		# between those with the prev gene on the same strand and those
		# with the previous gene on the opposite strand

		# If on the leading strand
		if cds_strand == '0':

			# If not the first cds
			if cds_id != 1:

				if cds_locus in filtered_genes[acc]:

					prev_cds_start = cds[cds_id-1][1]
					prev_cds_end = cds[cds_id-1][2]



					# If the cds start does not overlap the prev cds
					if int(cds_start) > int(prev_cds_end):


						distance_between = int(cds[cds_id][1]) - int(prev_cds_end)
						if distance_between:

							# Get the strand
							prev_cds_strand = cds[cds_id-1][3]

							total_with_cds_upstream += 1

							# If the same strand
							if prev_cds_strand == cds_strand:

								total_with_cds_upstream_on_same_strand += 1

								if cds_seq[3] == "A":
									cds_same.append(1)
									same_a_dist.append(distance_between)
									obs_cds_with_a_with_upstream_cds_same_strand += 1
									obs_cds_with_a_with_upstream_cds += 1
								else:
									cds_same.append(0)
									same_nota_dist.append(distance_between)
							else:
								if cds_seq[3] == "A":
									cds_opposite.append(1)
									opp_a_dist.append(distance_between)
									obs_cds_with_a_with_upstream_cds += 1
								else:
									cds_opposite.append(0)
									opp_nota_dist.append(distance_between)

		# If on the lagging strand
		elif cds_strand == "1":

			# If not the last cds
			if cds_id != len(cds):


				if cds_locus in filtered_genes[acc]:


					# Start and ends are opposite for genes on opposite strand

					prev_cds_start = cds[cds_id+1][1]
					prev_cds_end = cds[cds_id+1][2]
					prev_cds_strand = cds[cds_id+1][3]



					# If the cds start does not overlap the prev cds
					if cds_end < prev_cds_start:



						distance_between = int(prev_cds_start) - int(cds_end)
						if distance_between:

							total_with_cds_upstream += 1

							if prev_cds_strand == cds_strand:

								total_with_cds_upstream_on_same_strand += 1

								if cds_seq[3] == "A":
									cds_same.append(1)
									same_a_dist.append(distance_between)
									obs_cds_with_a_with_upstream_cds_same_strand += 1
									obs_cds_with_a_with_upstream_cds += 1
								else:
									cds_same.append(0)
									same_nota_dist.append(distance_between)
							else:
								if cds_seq[3] == "A":
									cds_opposite.append(1)
									opp_a_dist.append(distance_between)
									obs_cds_with_a_with_upstream_cds += 1
								else:
									cds_opposite.append(0)
									opp_nota_dist.append(distance_between)

	lead_dist_a = []
	lead_dist_nota = []
	mean_dist_lead = []
	count_a_lead = []
	lead_cds_len = 0

	mean_distance_same = []
	mean_distance_opposite = []

	dist_range_prop_lead = {}
	dist_range_prop_lag = {}

	for i in range(50,4550,50):
		dist_range_prop_lead[i] = []
		dist_range_prop_lag[i] = []

	for lead_id in lead_strand:

		cds_locus = lead_strand[lead_id][0]
		cds_start = lead_strand[lead_id][1]
		cds_end = lead_strand[lead_id][2]
		cds_strand = lead_strand[lead_id][3]
		cds_seq = lead_strand[lead_id][4]




		if lead_id != 1:
			if cds_locus in filtered_genes[acc]:

				prev_cds_end = lead_strand[lead_id-1][2]

				# If the cds does not overlap the upstream cds on the same strand
				if int(cds_start) > int(prev_cds_end):

					cds_distance = int(cds_start) - int(prev_cds_end)

					# Get the next multiple of 50 from the distance
					cds_dist_cat = ((cds_distance / 50) + 1) * 50


					mean_dist_lead.append(cds_distance)

					if cds_seq[3] == "A":
						lead_dist_a.append(cds_distance)
						count_a_lead.append(1)
						if cds_dist_cat in dist_range_prop_lead:
							dist_range_prop_lead[cds_dist_cat].append(1)

					else:
						lead_dist_nota.append(cds_distance)
						count_a_lead.append(0)
						if cds_dist_cat in dist_range_prop_lead:
							dist_range_prop_lead[cds_dist_cat].append(0)

	lag_dist_a = []
	lag_dist_nota = []
	mean_dist_lag = []
	count_a_lag = []
	lag_cds_len = 0


	for lag_id in lag_strand:

		cds_locus = lag_strand[lag_id][0]
		cds_start = lag_strand[lag_id][1]
		cds_end = lag_strand[lag_id][2]
		cds_strand = lag_strand[lag_id][3]
		cds_seq = lag_strand[lag_id][4]




		if lag_id != len(lag_strand):
			if cds_locus in filtered_genes[acc]:

				# Actually the cds end
				prev_cds_start = lag_strand[lag_id+1][2]

				# If the cds does not overlap the upstream cds on the same strand
				if int(cds_end) < int(prev_cds_start):

					cds_distance = int(prev_cds_start) - int(cds_end)

					cds_dist_cat = ((cds_distance / 50) + 1) * 50

					mean_dist_lag.append(cds_distance)


					if cds_seq[3] == "A":
						lag_dist_a.append(cds_distance)
						count_a_lag.append(1)
						if cds_dist_cat in dist_range_prop_lead:
							dist_range_prop_lag[cds_dist_cat].append(1)
					else:
						lag_dist_nota.append(cds_distance)
						count_a_lag.append(0)
						if cds_dist_cat in dist_range_prop_lead:
							dist_range_prop_lag[cds_dist_cat].append(0)




	prop_a_with_upstream = obs_cds_with_a_with_upstream_cds / float(total_with_cds_upstream)
	exp_a_upstream_same_strand = prop_a_with_upstream * total_with_cds_upstream_on_same_strand

	obs_exp = ((obs_cds_with_a_with_upstream_cds_same_strand - exp_a_upstream_same_strand)**2) / exp_a_upstream_same_strand
	chi_vals["all"].append(obs_exp)
	obs_exp_diff["all"].append(obs_cds_with_a_with_upstream_cds_same_strand - exp_a_upstream_same_strand)

	mean_dist_a_same = numpy.mean(same_a_dist)
	mean_dist_nota_same = numpy.mean(same_nota_dist)
	mean_dist_a_opp = numpy.mean(opp_a_dist)
	mean_dist_nota_opp = numpy.mean(opp_nota_dist)



	if len(cds_same) != 0:
		prop_a_same = sum(cds_same)/ float(len(cds_same))
	else:
		prop_a_same = 0

	if len(cds_opposite) != 0:
		prop_a_opposite = sum(cds_opposite) / float(len(cds_opposite))
	else:
		prop_a_opposite = 0
	output_file_line = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (acc,prop_a_same, prop_a_opposite,mean_dist_a_same,mean_dist_nota_same,mean_dist_a_opp,mean_dist_nota_opp,obs_cds_with_a_with_upstream_cds_same_strand,exp_a_upstream_same_strand,obs_exp)
	output_file.write(output_file_line)



	mean_dist_a_lead = numpy.mean(lead_dist_a)
	mean_dist_nota_lead = numpy.mean(lead_dist_nota)
	mean_dist_a_lag = numpy.mean(lag_dist_a)
	mean_dist_nota_lag = numpy.mean(lag_dist_nota)

	prop_a_lead = sum(count_a_lead)/float(len(count_a_lead))
	mean_dist_lead = numpy.mean(mean_dist_lead)
	prop_a_lag = sum(count_a_lag)/float(len(count_a_lag))
	mean_dist_lag = numpy.mean(mean_dist_lag)


	output_file2_line = "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (acc,mean_dist_a_lead,mean_dist_nota_lead,mean_dist_a_lag,mean_dist_nota_lag,prop_a_lead,mean_dist_lead,prop_a_lag,mean_dist_lag)
	output_file2.write(output_file2_line)


	for distance in sorted(dist_range_prop_lead.iterkeys()):
		if len(dist_range_prop_lead[distance]) != 0:
			prop_a_lead = sum(dist_range_prop_lead[distance])/float(len(dist_range_prop_lead[distance]))
		else:
			prop_a_lead = 0

		dist_a_prop_lead[distance].append(prop_a_lead)

	for distance in sorted(dist_range_prop_lag.iterkeys()):

		if len(dist_range_prop_lag[distance]) != 0:
			prop_a_lag = sum(dist_range_prop_lag[distance])/float(len(dist_range_prop_lag[distance]))
		else:
			prop_a_lag = 0

		dist_a_prop_lag[distance].append(prop_a_lag)


def upstream_analysis_distance(acc, cds, filtered_genes, output_file, output_file2, dist_a_prop_lead, dist_a_prop_lag, chi_vals, distance, obs_exp_diff):

	lead_strand = {}
	lag_strand = {}
	lead_count = 0
	lag_count = 0

	cds_same = []
	cds_opposite = []

	same_a_dist = []
	same_nota_dist = []
	opp_a_dist =[]
	opp_nota_dist = []


	obs_cds_with_a_with_upstream_cds_same_strand = 0
	total_with_cds_upstream = 0
	obs_cds_with_a_with_upstream_cds = 0
	total_with_cds_upstream_on_same_strand = 0

	for cds_id in cds:

		cds_locus = cds[cds_id][0]
		cds_start = cds[cds_id][1]
		cds_end = cds[cds_id][2]
		cds_strand = cds[cds_id][3]
		cds_seq = cds[cds_id][4]




		if cds_strand == "0":
			lead_count += 1
			lead_strand[lead_count] = [cds_locus,cds_start,cds_end,cds_strand,cds_seq]
		elif cds_strand == "1":
			lag_count += 1
			lag_strand[lag_count] = [cds_locus,cds_start,cds_end,cds_strand,cds_seq]


		# Look at the prop A of non-overlapping genes (either strand)
		# between those with the prev gene on the same strand and those
		# with the previous gene on the opposite strand

		# If on the leading strand
		if cds_strand == '0':

			# If not the first cds
			if cds_id != 1:

				if cds_locus in filtered_genes[acc]:

					prev_cds_start = cds[cds_id-1][1]
					prev_cds_end = cds[cds_id-1][2]



					# If the cds start does not overlap the prev cds
					if int(cds_start) > int(prev_cds_end):

						distance_between = int(cds[cds_id][1]) - int(prev_cds_end)

						if distance_between < distance:

							# Get the strand
							prev_cds_strand = cds[cds_id-1][3]

							total_with_cds_upstream += 1

							if cds_seq[3] == "A":
								obs_cds_with_a_with_upstream_cds += 1

							# If the same strand
							if prev_cds_strand == cds_strand:

								total_with_cds_upstream_on_same_strand += 1

								if cds_seq[3] == "A":
									cds_same.append(1)
									same_a_dist.append(distance_between)
									obs_cds_with_a_with_upstream_cds_same_strand += 1

								else:
									cds_same.append(0)
									same_nota_dist.append(distance_between)
							else:
								if cds_seq[3] == "A":
									cds_opposite.append(1)
									opp_a_dist.append(distance_between)

								else:
									cds_opposite.append(0)
									opp_nota_dist.append(distance_between)

		# If on the lagging strand
		elif cds_strand == "1":

			# If not the last cds
			if cds_id != len(cds):


				if cds_locus in filtered_genes[acc]:


					# Start and ends are opposite for genes on opposite strand

					prev_cds_start = cds[cds_id+1][1]
					prev_cds_end = cds[cds_id+1][2]
					prev_cds_strand = cds[cds_id+1][3]



					# If the cds start does not overlap the prev cds
					if cds_end < prev_cds_start:

						distance_between = int(prev_cds_start) - int(cds_end)

						if distance_between < distance:

							total_with_cds_upstream += 1

							if cds_seq[3] == "A":
									obs_cds_with_a_with_upstream_cds += 1

							if prev_cds_strand == cds_strand:

								total_with_cds_upstream_on_same_strand += 1

								if cds_seq[3] == "A":
									cds_same.append(1)
									same_a_dist.append(distance_between)
									obs_cds_with_a_with_upstream_cds_same_strand += 1
								else:
									cds_same.append(0)
									same_nota_dist.append(distance_between)
							else:
								if cds_seq[3] == "A":
									cds_opposite.append(1)
									opp_a_dist.append(distance_between)

								else:
									cds_opposite.append(0)
									opp_nota_dist.append(distance_between)

	lead_dist_a = []
	lead_dist_nota = []
	mean_dist_lead = []
	count_a_lead = []
	lead_cds_len = 0

	mean_distance_same = []
	mean_distance_opposite = []

	dist_range_prop_lead = {}
	dist_range_prop_lag = {}

	for i in range(50,4550,50):
		dist_range_prop_lead[i] = []
		dist_range_prop_lag[i] = []

	for lead_id in lead_strand:

		cds_locus = lead_strand[lead_id][0]
		cds_start = lead_strand[lead_id][1]
		cds_end = lead_strand[lead_id][2]
		cds_strand = lead_strand[lead_id][3]
		cds_seq = lead_strand[lead_id][4]




		if lead_id != 1:
			if cds_locus in filtered_genes[acc]:

				prev_cds_end = lead_strand[lead_id-1][2]

				# If the cds does not overlap the upstream cds on the same strand
				if int(cds_start) > int(prev_cds_end):

					cds_distance = int(cds_start) - int(prev_cds_end)

					# Get the next multiple of 50 from the distance
					cds_dist_cat = ((cds_distance / 50) + 1) * 50


					mean_dist_lead.append(cds_distance)

					if cds_seq[3] == "A":
						lead_dist_a.append(cds_distance)
						count_a_lead.append(1)
						if cds_dist_cat in dist_range_prop_lead:
							dist_range_prop_lead[cds_dist_cat].append(1)

					else:
						lead_dist_nota.append(cds_distance)
						count_a_lead.append(0)
						if cds_dist_cat in dist_range_prop_lead:
							dist_range_prop_lead[cds_dist_cat].append(0)

	lag_dist_a = []
	lag_dist_nota = []
	mean_dist_lag = []
	count_a_lag = []
	lag_cds_len = 0


	for lag_id in lag_strand:

		cds_locus = lag_strand[lag_id][0]
		cds_start = lag_strand[lag_id][1]
		cds_end = lag_strand[lag_id][2]
		cds_strand = lag_strand[lag_id][3]
		cds_seq = lag_strand[lag_id][4]




		if lag_id != len(lag_strand):
			if cds_locus in filtered_genes[acc]:

				# Actually the cds end
				prev_cds_start = lag_strand[lag_id+1][2]

				# If the cds does not overlap the upstream cds on the same strand
				if int(cds_end) < int(prev_cds_start):

					cds_distance = int(prev_cds_start) - int(cds_end)

					cds_dist_cat = ((cds_distance / 50) + 1) * 50

					mean_dist_lag.append(cds_distance)


					if cds_seq[3] == "A":
						lag_dist_a.append(cds_distance)
						count_a_lag.append(1)
						if cds_dist_cat in dist_range_prop_lead:
							dist_range_prop_lag[cds_dist_cat].append(1)
					else:
						lag_dist_nota.append(cds_distance)
						count_a_lag.append(0)
						if cds_dist_cat in dist_range_prop_lead:
							dist_range_prop_lag[cds_dist_cat].append(0)




	prop_a_with_upstream = obs_cds_with_a_with_upstream_cds / float(total_with_cds_upstream)
	exp_a_upstream_same_strand = prop_a_with_upstream * total_with_cds_upstream_on_same_strand


	obs_exp = ((obs_cds_with_a_with_upstream_cds_same_strand - exp_a_upstream_same_strand)**2) / exp_a_upstream_same_strand
	chi_vals[distance].append(obs_exp)
	obs_exp_diff[distance].append(obs_cds_with_a_with_upstream_cds_same_strand - exp_a_upstream_same_strand)

	if len(cds_same) != 0:
		prop_a_same = sum(cds_same)/ float(len(cds_same))
	else:
		prop_a_same = 0

	if len(cds_opposite) != 0:
		prop_a_opposite = sum(cds_opposite) / float(len(cds_opposite))
	else:
		prop_a_opposite = 0
	output_file_line = ",%s,%s,%s,%s,%s,%s" % (total_with_cds_upstream,obs_cds_with_a_with_upstream_cds,total_with_cds_upstream_on_same_strand,obs_cds_with_a_with_upstream_cds_same_strand,exp_a_upstream_same_strand,obs_exp)
	output_file.write(output_file_line)










def main():

	setupDirectory(output_directory)

	# Get a list of the filtered genes
	filtered_genes = get_filtered_genes(filtered_genes_file_path)

	# Get a list of table 4 genomes
	table4_genomes = get_table4_genomes(table4_genomes_path)

	output_file = open(output_file_path, "w")
	output_file.write("acc,prop_a_same,prop_a_opposite,mean_dist_a_same,mean_dist_nota_same,mean_dist_a_opp,mean_dist_nota_opp,obs_a_same,exp_a_same,((obs-exp)^2)/e")
	for distance in distances:
		output_line = ",cds_with_upstream_cds<%s,cds_with_a_with_upstream_cds<%s,total_cds_upstream_on_same_strand<%s,obs_cds_with_a_with_upstream_cds_on_same_strand<%s,exp_a_same<%s,((obs-exp)^2)/e<%s" % (distance,distance, distance, distance, distance, distance)
		output_file.write(output_line)
	output_file.write("\n")

	output_file2 = open(output2_file_path, "w")
	output_file2.write("acc,mean_dist_a_lead,mean_dist_nota_lead,mean_dist_a_lag,mean_dist_nota_lag,prop_a_lead,mean_dist_lead,prop_a_lag,mean_dist_lag\n")


	genome_count = 0
	dist_a_prop_lead = {}
	dist_a_prop_lag = {}

	obs_exp_diff = {}
	chi_vals = {}


	obs_exp_diff["all"] = []
	chi_vals["all"] = []

	for distance in distances:
		chi_vals[distance] = []
		obs_exp_diff[distance] = []


	for i in range(50,4550,50):
		dist_a_prop_lead[i] = []
		dist_a_prop_lag[i] = []

	for acc in filtered_genes:

		if acc not in table4_genomes:
		# if acc == "AE005174":
		# if acc in testFiles:
			genome_count+=1
			print "-" * 20
			print "Genome %d of %d\n%s\n" % (genome_count, len(filtered_genes), acc)

			# Extract the cds
			cds = read_genome(acc, filtered_genes, table4_genomes, genomes_cds_path)

			# Get codon frequency
			upstream_analysis_all(acc, cds, filtered_genes, output_file, output_file2, dist_a_prop_lead, dist_a_prop_lag, chi_vals, obs_exp_diff)
			for distance in distances:
				upstream_analysis_distance(acc, cds, filtered_genes, output_file, output_file2, dist_a_prop_lead, dist_a_prop_lag, chi_vals, distance, obs_exp_diff)

			output_file.write("\n")


	p_vals = {}
	obs_exp_means = {}
	mean_obs_exp = {}

	chi_sum_all = sum(chi_vals["all"])
	obs_exp_mean_all = sum(chi_vals["all"])/float(len(chi_vals["all"]))
	p_value_all = robjects.r('pchisq(%s,%s,lower.tail=FALSE)' % (chi_sum_all, len(chi_vals["all"])-1))[0]
	mean_obs_exp_all = sum(obs_exp_diff["all"])/float(len(obs_exp_diff["all"]))

	output_line = "chi_sum,,,,,,,,,%s" % (chi_sum_all)
	output_file.write(output_line)


	for distance in distances:

		chi_sum = sum(chi_vals[distance])
		obs_exp_mean = sum(chi_vals[distance])/float(len(chi_vals[distance]))
		p_value = robjects.r('pchisq(%s,%s,lower.tail=FALSE)' % (chi_sum, len(chi_vals[distance])-1))[0]
		mean_obs_exp[distance] = sum(obs_exp_diff[distance])/float(len(obs_exp_diff[distance]))


		p_vals[distance] = p_value
		obs_exp_means[distance] = obs_exp_mean

		output_line = ",,,,,,%s" % (chi_sum)
		output_file.write(output_line)
	output_file.write("\n")

	output_line = "p_val,,,,,,,,,%s" % (p_value_all)
	output_file.write(output_line)

	for distance in sorted(p_vals.iterkeys()):
		output_line = ",,,,,,%s" % (p_vals[distance])
		output_file.write(output_line)

	output_file.write("\n")
	output_line = "mean_((obs-exp)^2)/e,,,,,,,,,%s" % (obs_exp_mean_all)
	output_file.write(output_line)


	for distance in sorted(obs_exp_means.iterkeys()):
		output_line = ",,,,,,%s" % (obs_exp_means[distance])
		output_file.write(output_line)

	output_file.write("\n")
	output_line = "mean_obs-exp,,,,,,,,,%s" % (mean_obs_exp_all)
	output_file.write(output_line)


	for distance in sorted(mean_obs_exp.iterkeys()):
		output_line = ",,,,,,%s" % (mean_obs_exp[distance])
		output_file.write(output_line)




	output_file.close()
	output_file2.close()


if __name__ == "__main__":

	main()
