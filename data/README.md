## Output data
===

Output data for key results. Data.xlsx has compiled all key results to a single file. Descriptions of the corresponding raw file outputs can be seen below, with file descriptions an annotations found in data_file_descriptions.xlsx.

### Files
* data.xlsx

Contains grouped raw results. Also contains a summary page listing each sheet with the corresponding raw output output file generated by the scripts.

* data_file_descriptions.xlsx

Describes each of the raw output files found in data.xlsx. Column definitions are provided.

===
##### Raw output file descriptions

###### aod_scores
* _aod_all_aminos_all_genomes.csv: average of difference scores for each amino acid at the second position for all genomes
* _aod_second_amino.csv: average of difference scores for each amino acid when genomes are grouped by GC content

###### archaea
For each given site X in the CDSs:

* site_X_ratios_archaea.csv: enrichment ratios of site X, comparing use of nucleotide at that site with comparable sites throughout the transcriptome in archaea

* _chisqaure_site_4_archaea.csv: chi square test for A enrichment at the fourth site in archaea

###### expression
* _CAI_a_notA.csv: mean CAI values for CDSs with and without A at the fourth site
* _CAI_a_ratios.csv: genome mean CAI values with A4 ratios for comparison
* _high_low_CAI_compare.csv: genome proportion of CDSs with fourth site A dependent on CAI value
* start_codon_CAI_means.csv: mean CAI value for CDSs starting with ATG, GTG or TTG start codons for each genome

###### met_pair
* _met_codon_pair_summary.csv: summary of chi square values for each codon and each amino acid following an ATG in the CDSs
* _met_codon_pair.csv: observed, expected and chisq values for each codon following an ATG in the CDSs

###### nucleotide_conservation
* _codon_gc_variance.csv: variance in GC content at each position in the codons 2-30

###### ratio_testing
For each given site X in the CDSs:

* site_X_ratios.csv: enrichment ratios of site X, comparing use of nucleotide at that site with comparable sites throughout the transcriptome
* site_X_ratios_t4.csv: enrichment ratios of site X, comparing use of nucleotide at that site with comparable sites throughout the transcriptome in table 4 genomes

* _chisquare_site_4_no_overlap.csv: chi square test for A enrichment at the fourth site after accounting for four site overlaps
* _chisquare_site_4_t4.csv: chi square test for A enrichment at the fourth site in table 4 genomes
* _chisquare_site_4.csv: chi square test for A enrichment at the fourth site
* _site_4_no_overlap.csv: enrichment ratios of site 4, comparing use of nucleotide at that site with comparable sites throughout the transcriptome after removing four site overlaps
* _site_4_ratios_all_t4_genomes.csv: enrichment ratios of site 4, comparing use of nucleotide at that site with comparable sites throughout the transcriptome, for the expanded set of table 4 genomes


###### sd_sequences
* _sdSequenceAnalysis.csv: results comparing fourth site A use in the presence or absence of SD sequences, strength of SD sequences and distance away from the start codon of the SD sequence

###### second_amino
* _second_amino_ratios_t4.csv: ratio comparing second amino acid use with genome wide amino acid use for table 4 genomes for each amino acid
* _second_amino_ratios.csv: ratio comparing second amino acid use with genome wide amino acid use for each amino acid
* _second_amino_usage.csv: chisq calculations for each amino acid in the second position
* second_amino_usage_summary.csv: summary of the chisq calculations for each amino acid in the second position

###### species_compare
* _compare_ecoli_shigella.csv: proportion of CDSs with substitution in Shigella when compared with E. coli in the first position of codons 1-11 for each nucleotide

###### stop_codons
* stopDistanceGenome.csv: mean distance to the first +1 stop codon (not including those found site 2-4), mean distance to first +1 stop with fourth site A, mean distance to first +1 stop without fourth site A
* stopDistanceGenomeSecondStop.csv: mean distance to the second +1 stop codon (not including those found site 2-4), mean distance to second +1 stop with fourth site A, mean distance to second +1 stop without fourth site A
* stopDistanceGenomeThirdStop.csv: mean distance to the third +1 stop codon (not including those found site 2-4), mean distance to third +1 stop with fourth site A, mean distance to third +1 stop without fourth site A
* stopDistanceTable4Genome.csv: mean distance to the first +1 stop codon (not including those found site 2-4), mean distance to first +1 stop with fourth site A, mean distance to first +1 stop without fourth site A for table 4 genomes
* stopDistanceTable4GenomeSecondStop.csv: mean distance to the second +1 stop codon (not including those found site 2-4), mean distance to second +1 stop with fourth site A, mean distance to second +1 stop without fourth site A for table 4 genomes
* stopDistanceTable4GenomeThirdStop.csv: mean distance to the third +1 stop codon (not including those found site 2-4), mean distance to third +1 stop with fourth site A, mean distance to third +1 stop without fourth site A for table 4 genomes
