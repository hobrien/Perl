#!/usr/bin/perl -w

=head1 NAME

trunc.pl version 2, 22 April 2013

=head1 SYNOPSIS

cat XXX | trunc.pl start_coord end_coord 

=head1 DESCRIPTION

This will accept a sequence (or sequences) from STDIN, truncate it to the specified region, and
print the results.

-Can accept a plain sequence or fasta format.

-Negative numbers will be truncated from the end of the sequence.
-If the end of the region is not specified, it will truncate to the end of the sequence.

Options:

=head2 NOTES
This version had been rewritten to use the same basic code as rc.pl, which simplifies it 
considerably. Any edits to the main program here should be replicated in that script and vice versa.
=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use strict;
use warnings;
use Bio::SeqIO;

my $usage = "type perldoc trunc.pl for help";

my $seqin = Bio::SeqIO->new(
                            -file   => "cat |",
                            -format => 'fasta',
                            );

my $seqout = Bio::SeqIO->new(
                             -file => "| cat",
                             -format => 'fasta',
                             );

while ( my $seq = $seqin->next_seq ) {
  $seqout->write_seq($seq->translate);
}
