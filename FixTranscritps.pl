#!/usr/bin/perl -w

=head1 NAME

FixTranscripts.pl -v 2 23 April 2013

=head1 SYNOPSIS

cat gtf_file | FixTranscript.pl species > out_file 

=head1 DESCRIPTION
Modifies GTF file from Cufflinks:

-Adjacent exons on the same contig that have the same gene_name are modified to be 
different exons of the same transcript_id

-Non-adjacent exons that have the same gene_name are modified to have unique gene names
by appending a .# to the end of all but the first instance

-Exons without gene_names are named SPECIES_unique_#


=head2 NOTES

-The latter two changes are made so that htseq-count can be run with the option 
'-i gene_name', which will allow comparisons of homologs in different species.

-There are two problems with this approach, both of which will require better ortholog
identification:
  -genes that are shared among species but not present in the reference will not be compared
  -paralogs will be matched up arbitrarily for the cross-species comparisons

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################
use strict;
use warnings;
use FileParser qw(ParseGTF FlattenGTF);

my $species = shift;
unless ( $species) { die "No species name specified\n"; }

my $gene_id = 'init';
my $trans_id = 'init';
my $gene_name = 'init';
my $exon_num = 1;
my $contig = 'init';
my $tss_id = 'init';
my %gene_names;
my $unique_num = 0;

while (<>) {
  chomp;
  my %tags = % { ParseGTF($_) };
  unless ( $tags{'gene_id'} ) { die "missing tag gene_id\n"; }
  unless ( $tags{'transcript_id'} ) { die "missing tag transcript_id\n"; }
  unless ( $tags{'exon_number'} ) { die "missing tag exon_number\n"; }
  unless ( $tags{'gene_name'} ) {
    $unique_num ++;
    $tags{'gene_name'} = $species . '_unique_' . $unique_num;
  }
  if ( $tags{'gene_name'} eq $gene_name and $tags{'seqname'} eq $contig) {
    $exon_num ++;
    $tags{'exon_number'} = $exon_num;
    $tags{'gene_id'} = $gene_id;
    $tags{'transcript_id'} = $trans_id;
    $tags{'tss_id'} = $tss_id;
  }
  else {
    if ( $gene_names{$tags{'gene_name'}} ) {
      $tags{'gene_name'} = join(".", ($tags{'gene_name'}, $gene_names{$tags{'gene_name'}}));
      $gene_names{$tags{'gene_name'}} ++;
    }
    else {
      $gene_names{$tags{'gene_name'}} = 1;
    }
    $exon_num = 1;
    $gene_name = $tags{'gene_name'};
    $gene_id = $tags{'gene_id'};
    $trans_id = $tags{'transcript_id'};
    $contig = $tags{'seqname'};
    $tss_id = $tags{'tss_id'};
  }
  print FlattenGTF(\%tags), "\n";
}
