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
use File::Basename;
use Getopt::Long;
use FileParser qw(WriteFasta);

my $infilename;
my $outfilename;
my $min = 0;
my $max = 999999999999;

#Getopt::Long::Configure ("bundling");
GetOptions(
'infile:s' => \$infilename,
'outfile:s' => \$outfilename,
'minimum:s' =>\$min,
'maximum:s' => \$max
);

open(my $outfile, "> ".($outfilename || '-'))  or die;
open(my $infile, "< ".($infilename || '-'))  or die;

my $seq;
my $name;
my $gaps = 0;
while (<$infile>) {
  chomp;
  if ( $_ =~ /^>/ ) {
    if ( $seq and length($seq)-$gaps <= $max and length($seq)-$gaps >= $min ) { print $outfile "$name\n", WriteFasta($seq), "\n"; }
    $name = $_;
    $seq = '';
    $gaps = 0;
  }
  else {
    if ( $_ =~ /([^acgtnACGTN\-?])/ ) { warn "Non-standard character $1 in sequence $_\n"; }
    foreach ($_ =~ /-/g ) { $gaps ++; }
    $seq .= $_;
  }
}
if ( $name ) { print $outfile "$name\n"; }
if ( $seq and length($seq) - $gaps <= $max and length($seq) - $gaps >= $min ) { print $outfile WriteFasta($seq), "\n"; }

  
