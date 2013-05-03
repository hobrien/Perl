#!/usr/bin/perl -w

#a simple script to trim all numbers in a csv file to the specified number of digits
#only digits after the decimal are removed
#the zero before the decimal in number < 1 is not counted (and is added if not present)


use warnings;
use strict;

my $sd = shift;

while (<>) {
  chomp;
  my @fields;
  foreach ( split(/[\s]*,[\s]*/, $_) ) { push(@fields, SigDigits($_, $sd)); }
  print join(",", @fields), "\n";
}

sub SigDigits {
  my $num = shift;
  my $sd = shift;
  $sd ++;
  if ( $num =~ /^\d*\.?\d*$/ ) {
    $num =~ s/^0+//;
    if ( length($num) > $sd ) { 
      if ( $num < 10 ** $sd ) { $num = substr($num, 0, 5); }
      else { $num = int($num); }
    }
    $num =~ s/^\./0./;
    $num =~ s/\.$//;
  }
  return $num;
}