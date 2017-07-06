from Bio import SeqIO
import os

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



def run_test_genomes():

    # Open the output file and write the header line
    output_file = open('results_fourth_site.csv', 'w')
    header_line = 'file,seq_count,codon_count'
    for nt in ['A', 'C', 'G', 'T']:
        header_line += ',%s_fourth_count,%s_fourth_prop,%s_codon_count,%s_codon_prop,%s4_ratio' % (nt,nt,nt,nt,nt)
    header_line += '\n'
    output_file.write(header_line)



    # For each of the genomes in the directory containing the test genomes
    for genome in os.listdir('test_genomes/'):

        if genome != '.DS_Store':


            # Set counters to 0
            seq_count = 0
            codon_count =0
            fourth_count = {}
            codon_start_count = {}
            for nt in ['A', 'C', 'G', 'T']:
                fourth_count[nt] = 0
                codon_start_count[nt] = 0

            # Read the sequences using Biopython
            for seq_record in SeqIO.parse("test_genomes/" + genome, "fasta"):

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


            # Start the line to be output to file
            output_line = '%s,%s,%s' % (genome, seq_count, codon_count)


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

                # print('%s fourth count: %s' % (nt, fourth_count[nt]))
                # print('%s fourth prop: %s' % (nt, fourth_prop))
                # print('%s codon count: %s' % (nt, codon_start_count[nt]))
                # print('%s codon prop: %s' % (nt, codon_prop))
                # print('%s4 ratio: %s' % (nt, final_ratio))

            print('-' * 30)

            output_line += '\n'
            output_file.write(output_line)

    output_file.close()




def main():
    run_test_genomes()


#################################

if __name__ == "__main__":
	main()
