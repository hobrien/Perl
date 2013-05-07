#!/usr/bin/perl -w

=head1 NAME

RemoveDuplicates.pl version 1, 7 May 2013

=head1 SYNOPSIS

cat INFILE | RemoveDuplicates.pl > OUTFILE

or

RemoveDuplicates INFILE OUTFILE

=head1 DESCRIPTION

Deletes sequences with duplicate names

=head2 NOTES


=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use warnings;
use strict;
use Bio::SeqIO;

my $infilename = shift;
unless ($infilename ) { $infilename = "cat |";  }

my $outfilename = shift;
unless ($outfilename ) { $outfilename = "| cat";  }

my $infile = Bio::SeqIO->new('-file' => $infilename,
         '-format' => 'fasta') or die "could not open seq file $infilename\n";
my %seqs;	
while (my $seq = $infile->next_seq) {
  $seqs{$seq->display_id} = $seq;
}

my $outfile = Bio::SeqIO->new('-file' => $outfilename,
         '-format' => 'fasta') or die "could not open seq file $outfilename\n";

foreach (keys %seqs) {
  $outfile->write_seq($seqs{$_});
}