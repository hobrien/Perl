Perl
====
Perl scripts for genomic analysis:

AddMapMan.pl -Takes the output of MapMan annotation and adds results to a list of gene_IDs 

AlignHomologs.pl -Aligns, trims and makes trees for each cluster from GroupSequences.pl

AllSelaginella.pl -Runs a command on all KRAUS, MOEL, UNC and WILD datasets

CheckMono.py -Determine if Selaginella sequence tree is congruent with species tree

ClusterHomologs.pl -Clusters sequences by blast eValue using single-linkage clustering

CombineCounts.pl, -Combine the output from DESeq for each species, providing counts and significance

CompareSets.pl -Compares the values in the first column of files and produces Venn diagram code

ContigStats.pl -Prints basic info about the sequences in a file

ConvertSeq.pl - Convert among a variety of sequence and alignment file formats

Ensembl.pl -Testing the functionality of the Ensembl API

ExtractEntrez.pl -Search for genbank accession numbers and generat a entrez query to limit blast results

FilterSeq.pl -filtering sequences on min and/or max length

GetSeqs.pl -indexes fasta file and retrieves specified sequences

GetSeqsByGeneID.pl - Matches GeneIDs to contig names and outputs seqs matching list of GeneIDs

GroupSeqs.pl -Creates individual fasta files for each cluster identified by ClusterHomologs.pl

NonRef.pl - Creates a file sequences that do not match a list of names

rc.pl -Reverse-complement sequence

RemoveDuplicates.pl -Deletes sequences with duplicate names

ScafTranscritps.pl -Combines contigs that blast to the same reference gene

SigDigits.pl -a simple script to trim all numbers in a csv file to the specified number of digits

Temp - a folder for quick scripts that do not need to be added to the repo

TopHit.pl -Excludes sub-optimal hits for each query and sub-optimal queries for each hit

trunc.pl -truncates sequences from STDIN

=========================================================================================
Modules (Perl modules utilized by scripts in main folder):

FileParser.pl -subroutines to read and write GTF and (modified) blast tabular entries


=========================================================================================
t Folder of test datasets

test.bl - seqs 1,2,5 and 6 cluster (1 and 5 are opposite strand of 2 and 6) and 3 and 4 cluster (same strand)

test.fa - two sequences (each 1225bp plus 220 -'s, 26 n's), names correspond to test.bl

test.gtf -gtf with 5 exon features (3 match the same ref gene (2 consecutive), one does not match ref)

raw_seq.txt - one sequence (1225bp plus 220 -'s, 26 n's) with no fasta header


