#!/usr/bin/perl -w

=head1 NAME

AddMapMan.pl version 1, 17 April 2013

=head1 SYNOPSIS

AddMapMan.pl MapMan_file DEG_file > outfile

=head1 DESCRIPTION

Takes the output of MapMan annotation and adds results to a list of gene_IDs (such as 
from DEG analysis

=head2 NOTES

Currently searches for "XLOC_\d+" as the gene_id key (consistent with results from
TopHat). I will have to think about how to generalize the input gene_id.

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################
use warnings;
use strict;
use autodie qw(open close);

my %genes;

my $anno_filename = shift;
my $deg_filename = shift;

open(my $anno, $anno_filename);
while (<$anno>) {
  chomp;
  $_ =~ s/(\S)weakly similar to/$1; weakly similar to/g;
  $_ =~ s/(\S)very weakly similar to/$1; very weakly similar to/g;
  $_ =~ s/(\S)moderately similar to/$1; moderately similar to/g;
  $_ =~ s/(\S)highly similar to/$1; highly similar to/g;
  $_ =~ s/(\S)very highely similar to/$1; very highly similar to/g;
  $_ =~ s/\s+T$//;
  if ( $_ =~ s/(XLOC_\d+)\s+// ) {
    $genes{$1} = $_;
  }
}

open(my $deg, $deg_filename);
while ( <$deg> ) {
  chomp;
  print $_;
  if ( $_ =~ /(XLOC_\d+)/ and $genes{$1}) { print "\t", $genes{$1}; }
  print "\n";
}
  