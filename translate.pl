#!/usr/bin/perl -w

=head1 NAME

translate.pl version 3, 28 Jan 2014

=head1 SYNOPSIS

cat FILENAME | translate.pl

OR

translate.pl FILENAME

=head1 DESCRIPTION

This will accept a sequence (or sequences) from STDIN or one or more files, translate them
to amino acid sequences and print the results

-Can accept a plain sequence or fasta format.

=head2 NOTES
Switched back to manual parsing rather than use of SeqIO because I can't figure out how
to handle multiple formats (fasta/raw using SeqIO)

I tried to write a biopython version of this, but biopython complains if you try to 
translate a sequence with gaps, so this is more flexiable
=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use strict;
use warnings;
use Bio::Seq;

my $usage = "type perldoc trunc.pl for help";
my $name;
my $seq;
  
while (<>) {
  ($name, $seq) = parse_file($name, $seq, $_);
}
($name, $seq) = parse_file($name, $seq, '>');

sub parse_file {
  my $name = shift;
  my $seq = shift;
  my $line = shift;
  if ( $line =~ /^>/ ) {
    if ($name) {
      print $name;
    }
    if ($seq) {
      $seq =~ s/\s//g;
      print translate($seq)->seq, "\n";
    }
    $name = $line; 
    $seq = ''
  }
  else {  $seq .= $_; }
  return ($name, $seq)
}

sub translate {
  my $seq = shift;
  my $seq_obj = Bio::Seq->new(-seq => $seq,
                         -alphabet => 'dna' );
  return $seq_obj->translate();
}
