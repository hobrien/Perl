#!/usr/bin/perl -w

=head1 NAME

FilterSeq.pl -v 1 23 April 2013

=head1 SYNOPSIS

Filterseq.pl -i seq_file -o out_file --min min_size --max max_size
cat seq_file | Filterseq.pl --min min_size --max max_size > out_file 

=head1 DESCRIPTION
Filter sequences on length. Can specify min size and/or max size.
Can specify files for input and/or output or use STDIN/STDOUT

=head2 NOTES

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################


use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $infilename = "cat |";
my $outfilename = "| cat";
my $min = 0;
my $max = 999999999999;

#Getopt::Long::Configure ("bundling");
GetOptions(
'infile:s' => \$infilename,
'outfile:s' => \$outfilename,
'minimum:s' =>\$min,
'maximum:s' => \$max
);

unless ( $outfilename eq "| cat" ) { $outfilename = ">" . $outfilename; }

my $seqin = Bio::SeqIO->new(
                            -file   => $infilename,
                            -format => 'fasta',
                            );

my $seqout = Bio::SeqIO->new(
                             -file => $outfilename,
                             -format => 'Fasta',
                             );

while ( my $seq = $seqin->next_seq ) {
  my $gaps = 0;
  foreach ($seq->seq =~ /-/g ) { $gaps ++; };
  if ( $seq->length - $gaps <= $max and $seq->length - $gaps >= $min ) { $seqout->write_seq($seq); }
}
  
