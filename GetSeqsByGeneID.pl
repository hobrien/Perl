#!/usr/bin/perl -w

=head1 NAME

GetSeqsByGeneID.pl version 1, 17 April 2013

=head1 SYNOPSIS

GetSeqsByGeneID.pl MapMan_file DEG_file > outfile

=head1 DESCRIPTION

Takes GTF file and matches GeneIDs to contig names, then outputs sequences that correspond
to a list of GeneIDs.

=head2 NOTES

Skips first line. GeneID must be surrounded by quotes and in the first column.

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use warnings;
use strict;
use autodie qw(open close);
use FileParser qw(ParseGTF);
use Bio::Index::Fasta;
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";

#Prepare contig sequence index
my $contigfilename = shift;
my $contig_inx_name = $contigfilename . ".inx";
my $contig_inx = Bio::Index::Fasta->new(-filename => $contig_inx_name, -write_flag => 1);
$contig_inx->make_index($contigfilename);

my %contigs;
my $gtf_file = shift;
open(my $gtf, $gtf_file);
while (<$gtf>) {
  my %line = % { ParseGTF($_) };
  $contigs{$line{'gene_id'} } = $line{'seqname'};
}
 
while (<>) {
  if ( $. == 1 ) { next; }
  $_ =~ /"(\w+)"/;
  my $seq = $contig_inx->fetch($contigs{$1})  or die "could not find sequence $contigs{$1}\n";
  print ">", $1, "\n";
  print $seq->seq, "\n";
}  

