#!/opt/local/bin/perl -w

=head1 NAME

trunc.pl version 2, 22 April 2013

=head1 SYNOPSIS

cat XXX | trunc.pl start_coord end_coord 

=head1 DESCRIPTION

This will accept a sequence (or sequences) from STDIN, truncate it to the specified region, and
print the results.

-Can accept a plain sequence or fasta format.

- Sequences are bast 1

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

my $seq_start = shift or die $usage;
my $seq_end = shift;
$seq_start -= 1;

my $seq;
while (<>) {
  chomp;
  if ( $_ =~ /^>/ ) {
    if ( $seq ) { print substr($seq, $seq_start, $seq_end), "\n"; }
    print "$_\n";
    $seq = '';
  }
  else {
    if ( $_ =~ /([^acgtnACGTN])/ ) { warn "Non-standard character $1 in sequence\n"; }
    $seq .= $_;
  }
}
print substr($seq, $seq_start, $seq_end), "\n";

