#!/usr/bin/perl -w

=head1 NAME

GroupSeqs.pl version 1, 7 May 2013

=head1 SYNOPSIS

cat Cluster_IDs | GroupSeqs.pl

=head1 DESCRIPTION

Creates individual fasta files for each cluster identified by ClusterHomologs.pl


=head2 NOTES

VERY IMPORTANT to run this script from within a new folder because it generates 
thousands of outfiles

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Index::Fasta;
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";


my $infilename = shift;
my $inx_name = $infilename . ".inx";
my $inx = Bio::Index::Fasta->new(-filename => $inx_name, -write_flag => 1);

while (<>) {
  chomp;
  my @fields = split(/\s*\t\s*/, $_);
  my $seq = $inx->fetch($fields[0]) or die "could not find sequence $fields[0]\n";
  $seq = $seq->trunc($fields[3], $fields[4]);
  my $outfilename = $fields[1] . '.fa';
  my $outfile = Bio::SeqIO->new('-file' => ">>$outfilename",
         '-format' => 'fasta') or die "could not open seq file $outfilename\n";
  if ( $fields[2] == 1 ) { $outfile->write_seq($seq); }
  else { $outfile->write_seq($seq->revcom); }
}