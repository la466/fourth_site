#!/usr/bin/python

# Script number: 			4.6
# File: 					6 of 6
# Prerequisite script(s):	1.1, 1.2, 2.1, 3.1, 3.2
# Prerequisite file(s):		_filteredGenes.csv
# Description: 				Fourth site use for all table 4 genomes
# Output file(s):			site_4_ratios_all_t4_genomes.csv


from Bio import SeqIO
import os
import multiprocessing as mp


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


def check_seq(seq):

    # Set the seq_pass to true
    seq_pass = True

    # Check the sequence has length multiple of 3
    if len(seq) % 3 != 0:
        seq_pass = False

    # Check the sequence starts with an NTG start codon
    if seq[:3] not in ['ATG', 'CTG', 'GTG', 'TTG']:
        seq_pass = False

    # Check the sequence has a table 4 stop codon
    if seq[-3:] not in ['TAA', 'TAG']:
        seq_pass = False

    # Check there are no in frame stops
    for i in range(0, len(seq)-3, 3):
        if seq[i:i+3] in ['TAA', 'TAG']:
            seq_pass = False

    # Check that the sequence only contains A, C, T or G
    for i in range(0,len(seq)):
        if seq[i] not in ['A', 'C', 'T', 'G']:
            seq_pass = False

    # Return the filter check
    return(seq_pass)



def run_genomes():

    # Open the output file and write the header line
    output_file = open('outputs/ratio_testing/site_4_ratios_all_t4_genomes.csv', 'w')
    header_line = 'acc,genus,gc,seq_count,codon_count'
    for nt in ['A', 'C', 'G', 'T']:
        header_line += ',%s_fourth_count,%s_fourth_prop,%s_codon_count,%s_codon_prop,%s4_ratio' % (nt,nt,nt,nt,nt)
    header_line += '\n'
    output_file.write(header_line)

    genomes = []

    genome_count = 0
    # For each of the genomes in the directory containing the test genomes
    for genome in os.listdir('genome_extractions/t4/'):
        genome_count += 1
        if genome_count <= 20:
            genomes.append(genome)




    results = run_in_parralell(genomes, [], calc_ratios)
    for result in results:
        outputs = result.get()
        for output in outputs:
            output_file.write(output)

    output_file.close()

def calc_ratios(genomes):

    outputs = []

    for genome in genomes:

        # Set counters to 0
        seq_count = 0
        codon_count =0
        fourth_count = {}
        codon_start_count = {}
        nt_count = 0
        gc_nt_count = {}

        got_genus = False

        for nt in ['A', 'C', 'G', 'T']:
            fourth_count[nt] = 0
            codon_start_count[nt] = 0
            gc_nt_count[nt] = 0

        # Read the sequences using Biopython
        for seq_record in SeqIO.parse("genome_extractions/t4/" + genome, "fasta"):

            if not got_genus:
                desc = seq_record.description
                splits = desc.split(';')
                genus = splits[0]
                got_genus = True


            # Return the CDS
            seq = seq_record.seq

            # Check the sequences to pass the filtering criteria
            seq_pass = check_seq(seq)

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


                for i in range(0, len(seq)-3):
                    gc_nt_count[seq[i]] += 1
                nt_count += len(seq[:-3])

        gc = (gc_nt_count['G'] + gc_nt_count['C'])/nt_count

        # Start the line to be output to file
        output_line = '%s,%s,%s,%s,%s' % (genome.strip('.txt'),genus,gc, seq_count, codon_count)


        print(genome)

        for nt in sorted(fourth_count):

            # Calculate the proportion of CDSs using nt at fourth site
            fourth_prop = fourth_count[nt] / seq_count

            # Calcualte the proportion of genome codons start with nt
            codon_prop = codon_start_count[nt] / codon_count

            # Calculate N4 ratio
            final_ratio = fourth_prop / codon_prop

            # Add nt information to output line
            output_line += ',%s,%s,%s,%s,%s' % (fourth_count[nt], fourth_prop, codon_start_count[nt], codon_prop, final_ratio)



        print('-' * 30)



        output_line += '\n'
        outputs.append(output_line)
        # output_file.write(output_line)

    return(outputs)




def main():
    run_genomes()


#################################

if __name__ == "__main__":
	main()
