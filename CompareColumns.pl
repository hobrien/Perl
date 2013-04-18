#!/usr/bin/perl -w

=head1 NAME

CompareColumns.pl version 1, 17 April 2013

=head1 SYNOPSIS

CompareColumns.pl file1 index1 file2 index2 file3 index3

=head1 DESCRIPTION

Compares the values in the columns specified (index) and produces prints code that can be 
used to generate a Venn diagram in R depicting the overlaps.

NB: uses 1-based indices

Can be extended to allow variable numbers of comparisons (currently set to 3) and to
automatically generate the figure in R

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-utoronto-dot-caE<gt>

=cut


####################################################################################################
#see http://statistics.ats.ucla.edu/stat/r/faq/venn.htm for instructions on preparing the input file
#use library(VennDiagram) in R
use strict;
use warnings;


my $file1name = shift;
my $index1 = shift;
$index1 --;
my $file2name = shift;
my $index2 = shift;
$index2 --;
my $file3name = shift;
my $index3 = shift;
$index3 --;

my %file1_genes;
my $file1unique = 0;
open(my $file1, $file1name);
while (<$file1>) {
  chomp;
  if ( $. ==0 or $_ =~ /^\s*$/ or $_ =~ /Not\s+in\s+POOL/) { print "skipping $_\n"; next; }
  my @fields = split(/\s*\t\s*/, uc($_));
  if ( $fields[$index1] =~/^none$/i ) { $file1unique ++; next; }
  if ( $fields[$index1] =~/^na$/i ) { $file1unique ++; next; }
  $file1_genes{$fields[$index1]} = 1;
}
close($file1);

my %file2_genes;
my $file2unique = 0;
open(my $file2, $file2name);
while (<$file2>) {
  chomp;
  if ( $. ==0 or $_ =~ /^\s*$/ or $_ =~ /Not\s+in\s+POOL/) { print "skipping $_\n"; next; }
  my @fields = split(/\s*\t\s*/, uc($_));
  if ( $fields[$index2] =~/^none$/i ) { $file2unique ++; next; }
  if ( $fields[$index2] =~/^na$/i ) { $file2unique ++; next; }
  $file2_genes{$fields[$index2]} = 1;
}
close($file2);

my %file3_genes;
my $file3unique = 0;
open(my $file3, $file3name);
while (<$file3>) {
  chomp;
  if ( $. ==0 or $_ =~ /^\s*$/ or $_ =~ /Not\s+in\s+POOL/) { print "skipping $_\n"; next; }
  my @fields = split(/\s*\t\s*/, uc($_));
  if ( $fields[$index3] =~/^none$/i ) { $file3unique ++; next; }
  if ( $fields[$index3] =~/^na$/i ) { $file3unique ++; next; }
  $file3_genes{$fields[$index3]} = 1;
}
close($file3);

my $file123 =0;
my $file12 = 0;
my $file13 = 0;
my $file23 = 0;
my $area1 = 0;
my $area2 = 0;
my $area3 = 0;

foreach( keys %file1_genes ) {
  if ( $file2_genes{$_} ) {
    if ( $file3_genes{$_} ) {
      $file123 ++;
      $file12 ++;
      $file13 ++;
      $file23 ++;
      $area1 ++;
      $area2 ++;
      $area3 ++;
      delete $file3_genes{$_};
    }
    else {
      $file12 ++;
      $area1 ++;
      $area2 ++;
    }
    delete $file2_genes{$_};
  }
  elsif ( $file3_genes{$_} ) {
    $file13 ++;
      $area1 ++;
      $area3 ++;
    delete $file3_genes{$_};
  }
  else {
    $file1unique ++;
  }
  delete $file1_genes{$_};
}

foreach( keys %file2_genes ) {
  if ( $file3_genes{$_} ) {
    $file23 ++;
      $area2 ++;
      $area3 ++;
    delete $file3_genes{$_};
  }
  else {
    $file2unique ++;
  }
  delete $file2_genes{$_};
}
$file3unique += scalar( keys %file3_genes );

$area1 += $file1unique;
$area2 += $file2unique;
$area3 += $file3unique;

$file1name =~ s/.txt//;
$file2name =~ s/.txt//;
$file3name =~ s/.txt//;

print "ven.plot<-draw.triple.venn(
area1 = $area1,
area2 = $area2,
area3 = $area3,
n12 = $file12,
n23 = $file23,
n13 = $file13,
n123 = $file123,
category = c('$file1name', '$file2name', '$file3name'),
fill = c('blue', 'red', 'green'),
lty = 'blank',
cex = 2,
cat.cex = 1.75,
cat.col = c('blue', 'red', 'green'),
margin=0.05,
fontfamily='sans',
cat.fontfamily='sans'
)
";
