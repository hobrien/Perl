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
use Bio::Seq;

my $usage = "type perldoc trunc.pl for help";

my $name = 1;
if ( $ARGV[0] eq "-s" ) { $name = 0; shift; }
  
my $seq_start = shift or die $usage;
my $seq_end;
unless ( $seq_end = shift ) { $seq_end = 'end'; }
my $seq;
while (<>) {
  chomp;
  if ( $_ =~ /^>/ ) {
    if ( $seq ) { print TruncSeq($seq, $seq_start, $seq_end), "\n"; }
    print "$_\n";
    $seq = '';
  }
  else {
    if ( $_ =~ /([^acgtnACGTN])/ ) { warn "Non-standard character $1 in sequence $_\n"; }
    $seq .= $_;
  }
}
print TruncSeq($seq, $seq_start, $seq_end), "\n";


sub TruncSeq {
  my $seq = shift;
  my $seq_start = shift;
  my $seq_end = shift;
  if ( $seq_start < 0 ) { $seq_start = length($seq) - $seq_start + 1; }
  if ( $seq_end eq 'end' ) { $seq_end = length($seq); }
  elsif ( $seq_end < 0 ) { $seq_end = length($seq) - $seq_end + 1; }
  if ( $seq_end < $seq_start ) { my $temp = $seq_start; $seq_start = $seq_end; $seq_end = $temp; }
  unless ($seq_start == 0 ) { $seq_start --; }
#  print "start: $seq_start, length: ", $seq_end, " - ", $seq_start, " = ", $seq_end - $seq_start, "\n";
  return substr($seq, $seq_start, $seq_end - $seq_start);
}
