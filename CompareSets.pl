#!/usr/bin/perl -w

=head1 NAME

CompareSets.pl version 1, 24 June 2012

=head1 SYNOPSIS

CompareSets.pl file1 file2 [file3 file4] 

=head1 DESCRIPTION

Compares the values in the first column of three or four datasets and prints an R
command to generate a Venn Diagram of overlaps.

=head2 NOTES

use library(VennDiagram) in R to generate the plot

Options:

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use strict;
use warnings;
use List::Compare;
use List::MoreUtils qw(uniq);
use File::Basename;

my @sets;
my @files; 
foreach(@ARGV) {
  #$_ =~ /(\w+)/;
  (my $name, my $path, my $in_ext) = fileparse($_, qr/\.[^.]*/);
  push(@files, $name);
  open(my $in, "<", $_);
  my @names;
  while(<$in>) {
    my @fields = split(/\s+/, lc($_));
    if ( @fields == 0 ) { next; }
    if ( $fields[0] =~ /\S/ ) {
      $fields[0] =~ s/^\s+//;     #strip out leading whitespace
      $fields[0] =~ s/\s+$//;     #strip out leading whitespace      
      push(@names, $fields[0]); }
  }
  @names = uniq(@names);
  push(@sets, \@names);
}
print "# drawing Venn Diagram for total of ", scalar(List::Compare->new(@sets)->get_union), " items (type library(VennDiagram) to view in R)\n";

if ( @sets == 2 ) {
  print "ven.plot<-draw.pairwise.venn(\n";
}
elsif ( @sets == 3 ) {
  print "ven.plot<-draw.triple.venn(\n";
}
elsif ( @sets == 4 ) {
  print "venn.plot <- draw.quad.venn(\n";
}
elsif ( @sets == 5 ) {
  print "venn.plot <- draw.quintuple.venn(\n";
}
else { die "This script currently only produces pariwise to quintuple Venn plots\n"; }

for (my $i = 0; $i < @sets; $i ++) { 
  print "  area", $i+1, " = ", scalar(@{$sets[$i]}), ",\n";
}

for (my $i = 0; $i < @sets - 1; $i ++) {
  for (my $j = $i + 1; $j < @sets; $j ++ ) {
    my $lc= List::Compare->new($sets[$i], $sets[$j]);
    if ( @sets == 2 ) {     
      print "  cross.area = ", scalar($lc->get_intersection), ",\n";
    }
    else {
      print "  n", $i+1, $j+1, " = ", scalar($lc->get_intersection), ",\n";
    }
  }
}

for (my $i = 0; $i < @sets - 2; $i ++) {
  for (my $j = $i + 1; $j < @sets - 1; $j ++ ) {
    for ( my $k = $j + 1; $k < @sets; $k ++ ) {
      my $lc= List::Compare->new($sets[$i], $sets[$j], $sets[$k]);
      print "  n", $i+1, $j+1, $k+1, " = ", scalar($lc->get_intersection), ",\n";
    }
  }
}

for (my $i = 0; $i < @sets - 3; $i ++) {
  for (my $j = $i + 1; $j < @sets - 2; $j ++ ) {
    for ( my $k = $j + 1; $k < @sets - 1; $k ++ ) {
      for ( my $l = $k + 1; $l < @sets; $l ++ ) {      
        my $lc= List::Compare->new($sets[$i], $sets[$j], $sets[$k], $sets[$l]);
        print "  n", $i+1, $j+1, $k+1, $l+1, " = ", scalar($lc->get_intersection), ",\n";
      }
    }
  }
}

for (my $i = 0; $i < @sets - 4; $i ++) {
  for (my $j = $i + 1; $j < @sets - 3; $j ++ ) {
    for ( my $k = $j + 1; $k < @sets - 2; $k ++ ) {
      for ( my $l = $k + 1; $l < @sets -1; $l ++ ) {
        for ( my $m = $l + 1; $m < @sets; $m ++ ) {    
          my $lc= List::Compare->new($sets[$i], $sets[$j], $sets[$k], $sets[$l], $sets[$m]);
          print "  n", $i+1, $j+1, $k+1, $l+1, $m+1, " = ", scalar($lc->get_intersection), ",\n";
        }
      }
    }
  }
}

print "  category = c('", join("','", @files), "'),\n";
if ( @sets == 2) {
  print "  fill = c('blue', 'red'),
  cat.col = c('blue', 'red'),\n";
}  
elsif ( @sets == 3 ) {
  print "  fill = c('blue', 'red', 'green'),
  cat.col = c('blue', 'red', 'green'),\n";
}
elsif ( @sets == 4 ) {
  print "  fill = c('orange', 'red', 'green', 'blue'),
  cat.col = c('orange', 'red', 'green', 'blue'),\n";
}
else {
  print "  fill = c('orange', 'red', 'green', 'blue', 'purple'),
  cat.col = c('orange', 'red', 'green', 'blue', 'purple'),\n";
}

print "  lty = 'blank',
  cex = 2,
  cat.cex = 1.75,
  margin=0.05,
  fontfamily='sans',
  cat.fontfamily='sans'
 )\n";

	
