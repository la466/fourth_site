#!/usr/bin/python

# Script number: 			18.2
# File: 					2 of 2
# Prerequisite script(s):
# Prerequisite file(s):
# Description: 				Get the ratios for protists
# Output file(s):			_protists.csv

from Bio import SeqIO
import os
import multiprocessing as mp
import numpy
import re
import scipy.stats


trans_table = {}
trans_table[4] = ['TAA', 'TAG']
trans_table[6] = ['TGA']
trans_table[11] = ['TAA', 'TAG', 'TGA']

table6 = ['oxytricha', 'stylonychia', 'paramecium', 'tetrahymena', 'oxytrichidae']



def create_directory(direc):

    if not os.path.exists:
        os.makedirs(direct)


def run_in_parralell(input_list, args, function_to_run, workers = None, onebyone = False):

    if not workers:
        workers = (mp.cpu_count()) - 1

    if not onebyone:
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        chunk_list = input_list

    pool = mp.Pool(workers)
    results = []

    for i in chunk_list:
        current_args = args.copy()
        new_args = [i]



        for arg in current_args:
            # print(arg)
            new_args.append(arg)

        process = pool.apply_async(function_to_run, new_args)
        results.append(process)

    pool.close()
    pool.join()

    return(results)


def get_stops(organism):

    genus = re.findall(r'^[^_]+(?=_)', organism)[0]
    if genus.lower() in table6:
        table = 6
    else:
        table = 11

    return(trans_table[table])

def check_seq(seq, stops):

    # Set the seq_pass to true
    seq_pass = True

    error_codes = []

    # Check the sequence has length multiple of 3
    if len(seq) % 3 != 0:
        seq_pass = False
        error_codes.append(1)

    # Check the sequence starts with an NTG start codon
    if seq[:3] not in ['ATG', 'CTG', 'GTG', 'TTG']:
        seq_pass = False
        error_codes.append(2)

    # Check the sequence has a table 4 stop codon
    if seq[-3:] not in stops:
        seq_pass = False
        error_codes.append(3)

    # Check there are no in frame stops
    for i in range(0, len(seq)-3, 3):
        if seq[i:i+3] in stops:
            seq_pass = False
            error_codes.append(4)

    # Check that the sequence only contains A, C, T or G
    for i in range(0,len(seq)):
        if seq[i] not in ['A', 'C', 'T', 'G']:
            seq_pass = False
            error_codes.append(5)

    # Return the filter check
    return(seq_pass, error_codes)



def run_genome(genomes):

    outputs = []

    for genome in genomes:

        print(genome)

        stops = get_stops(genome)



        # Set counters to 0
        seq_count = 0
        codon_count =0
        fourth_count = {}
        codon_start_count = {}
        for nt in ['A', 'C', 'G', 'T']:
            fourth_count[nt] = 0
            codon_start_count[nt] = 0


        # Read the sequences using Biopython
        for seq_record in SeqIO.parse("protists/" + genome, "fasta"):




            # Return the CDS
            seq = seq_record.seq

            # Check the sequences to pass the filtering criteria
            seq_pass, error_codes = check_seq(seq, stops)

            # print(error_codes)
            # print(seq)

            # If the sequence passed the filtering
            if seq_pass:

                # Increment the sequence count
                seq_count += 1


                # Increment the fourth site count
                fourth_count[seq[3]] += 1

                # For each codon in the sequence (excluding the stop codon)
                for i in range(0, len(seq)-3, 3):
                    codon = seq[i:i+3]

                    # Increment the codon count
                    codon_count += 1

                    # Increment the codon start count
                    codon_start_count[codon[0]] += 1

        output = [genome, seq_count, codon_count, fourth_count, codon_start_count]
        outputs.append(output)



    return(outputs)



def divide(num,den):

    if (den) != 0:
        result = num/float(den)
    else:
        result = 0

    return(result)



def get_single_genus(genomes):

    genuses = {}
    single_genuses = []

    for genome in genomes:
        genus = re.findall(r'^[^_]+(?=_)', genome)[0]
        if genus.lower() not in genuses:
            genuses[genus.lower()] = [genome]
        else:
            genuses[genus.lower()].append(genome)

    for genus in genuses:
        numpy.random.seed()
        choice = numpy.random.choice(genuses[genus])
        single_genuses.append(choice)


    return(single_genuses)


def run():

    create_directory('outputs/')
    create_directory('outputs/protists')


    genomes = []
    genome_count = 0

    # For each of the genomes in the directory containing the test genomes
    for genome in os.listdir('protists/'):
        genome_count += 1
        if genome_count:
            genomes.append(genome)

    single_genus = get_single_genus(genomes)

    # # results = run_genome(genomes)
    results = run_in_parralell(single_genus, [], run_genome)


    # Open the output file and write the header line
    output_file = open('outputs/protists/protists.csv', 'w')
    header_line = 'file,seq_count,codon_count'
    for nt in ['A', 'C', 'G', 'T']:
        header_line += ',%s_fourth_count,%s_fourth_prop,%s_codon_count,%s_codon_prop,%s4_ratio' % (nt,nt,nt,nt,nt)
    header_line += '\n'
    output_file.write(header_line)


    print('generating results')
    for result in results:

        outputs = result.get()
        for output in outputs:

            genome = output[0]
            seq_count = output[1]
            codon_count = output[2]
            fourth_count = output[3]
            codon_start_count = output[4]


            # Start the line to be output to file
            output_line = '%s,%s,%s' % (genome, seq_count, codon_count)


            for nt in sorted(fourth_count):

                # Calculate the proportion of CDSs using nt at fourth site
                fourth_prop = divide(fourth_count[nt], seq_count)

                # Calcualte the proportion of genome codons start with nt
                codon_prop = divide(codon_start_count[nt],codon_count)

                # Calculate N4 ratio
                final_ratio = divide(fourth_prop, codon_prop)





                # Add nt information to output line
                output_line += ',%s,%s,%s,%s,%s' % (fourth_count[nt], fourth_prop, codon_start_count[nt], codon_prop, final_ratio)





            output_line += '\n'

            output_file.write(output_line)


    output_file.close()




def main():
    run()


#################################

if __name__ == "__main__":
	main()
