Perl
====
Perl scripts for genomic analysis:

AddMapMan.pl -Takes the output of MapMan annotation and adds results to a list of gene_IDs 

CompareColumns.pl -Compares the values in the specified columns and produces Venn diagram code

ContigStats.pl -Prints basic info about the sequences in a file.

Ensembl.pl -Testing the functionality of the Ensembl API.

GetSeqsByGeneID.pl - Matches GeneIDs to contig names and outputs seqs matching list of GeneIDs

NonRef.pl - Creates a file sequences that do not match a list of names

ScafTranscritps.pl -Combines contigs that blast to the same reference gene.

TopHit.pl -Excludes sub-optimal hits for each query and sub-optimal queries for each hit

=========================================================================================
Modules (Perl modules utilized by scripts in main folder):

FileParser.pl -subroutines to read and write GTF and (modified) blast tabular entries