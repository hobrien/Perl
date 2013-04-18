#!/usr/bin/perl -w

=head1 NAME

TopHit.pl version 2, 18 April 2013

=head1 SYNOPSIS

cat BLASTRESULT | TopHit.pl > OUTFILE

=head1 DESCRIPTION

Excludes all but the top hit for each query. In cases where 2 or more queries have the
same top hit, it also exlcudes all but the one with the highest score.

Result is one blast result for each query and for each hit.

=head2 NOTES

Does not combine multiple hits to the same hit sequence, which will definitley be 
necessary when applying it to scaffolded results. The code to deal with this is in
CombineHits.pl, which can be run on blast results prior to running this script, or the
code can be added here.

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use warnings;
use strict;
use Bio::Index::Fasta;
use Bio::SeqIO;
use FileParser qw(ParseBlast FlattenBlast);
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";

my $usage = "cat BLASTFILE | TopHit.pl > OUTFILE";

#parse blast hits and make a hash of tops hits
my $query_name = 'Init';
my %top_hits;

while (<>) {
  chomp;
  $_ =~ s/,\s*/\t/g;
  if ( $_ =~ /^#/ or $_ =~ /ADH\d{5}\t/) { next; }
  my %result = % {ParseBlast($_) };
  if ( $result{'query_name'} eq $query_name ) { next; } #not the top hit for a contig
  $query_name = $result{'query_name'};
  if ( $top_hits{$result{'hit_name'}} ) {               #other query matches same hit
    my %previous_result = % { ParseBlast($top_hits{$result{'hit_name'}}) };
    unless ( $result{'score'} > $previous_result{'score'} ) { next; }
  }
  $top_hits{$result{'hit_name'}} = FlattenBlast(\%result);
}  

print "Query_ID	Query_length	Hit_ID	Hit_length	Percent_ID	Aligned_length	Mismatch_count	Gap_count	Query_start	Query_end	Query_frame	Hit_start	Hit_end	Hit_frame	Evalue	Score\n";
foreach ( keys %top_hits ) {
  print $top_hits{$_}, "\n";
}