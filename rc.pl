#!/usr/bin/perl -w

=head1 NAME

rc.pl version 1, 9 June 2010

=head1 SYNOPSIS

cat XXX | rc.pl

=head1 DESCRIPTION

This will accept a sequence (or sequences) from STDIN, reverse complement it, and write the results
to STDOUT

can be plain sequence or fasta format.
=head2 NOTES

Use RevCom for sequences in files

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-utoronto-dot-caE<gt>

=cut
####################################################################################################

use strict;
use warnings;


my $usage = "Type perldoc rc.pl for help\n";

my $seq;
while (<>) {
  chomp;
  if ( $_ =~ /^>/ ) {
    if ( $seq ) { print RevCom($seq), "\n"; }
    print "$_\n";
    $seq = '';
  }
  else {
    if ( $_ =~ /([^acgtnACGTN])/ ) { warn "Non-standard character $1 in sequence\n"; }
    $seq .= $_;
  }
}
print RevCom($seq), "\n";

sub RevCom {
  my $seq = shift;
  $seq =~  tr/ACGT/TGCA/;
  $seq = reverse($seq);
  return $seq;
}
