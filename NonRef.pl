#!/usr/bin/perl -w

=head1 NAME

NonRef.pl version 1, 8 April 2013

=head1 SYNOPSIS

cut -f1 scaf_log_file | NonRef.pl seq_file

=head1 DESCRIPTION

Produces a file of sequences that do not match the reference. 

Designed to work with the output of ScafTranscripts.pl, but any other list of sequence
names to be excluded could be piped in 


=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut


####################################################################################################

use warnings;
use strict;
use Bio::SeqIO;

my %ref_genes;
my $infilename = shift;

while (<>) {
  chomp;
  $ref_genes{$_} = 1;
}
my $infile = Bio::SeqIO->new('-file' => $infilename,
         '-format' => 'fasta') or die "could not open seq file $infilename\n";

$infilename =~ s/\.fa/_nonref.fa/;
my $outfile = Bio::SeqIO->new('-file' => ">$infilename",
         '-format' => 'fasta') or die "could not open seq file $infilename\n";

while ( my $seq = $infile->next_seq ){
  unless ( $ref_genes{$seq->display_id} or $seq->length < 500 ) { $outfile->write_seq($seq); }
}
