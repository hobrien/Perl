#!/usr/bin/perl -w

#This is being rewritten to make a single exon feature for the length of the contig with a
#single cds feature the length of the top blast hit. Chimeric transcripts are going to be
#a big problem for this approach, but I can check for them independently (by comparing cds
#and exon length) and break them up if necessary before doing the expression analysis

use warnings;
use strict;
use FileParser qw(ParseBlast FlattenGTF);

my $previous_query = 'init';
my %exon = (
  'source' => 'Trinity',
  'feature' => 'exon',
  'frame' => '.',
  'score' => '.'
);
my %cds = (
  'source' => 'JGI',
  'feature' => 'CDS',
  'frame' => 1
);
while (<>) {
  chomp;
  my %blast = % { ParseBlast($_) };
  if ( $blast{'percent'} < 70 ) { #skip low-similarity hits
    next;
  }
    #skip secondary hits from same query
  if ( $blast{'query_name'} eq $previous_query ) { #skip secondary hits
    next; 
  }
  $previous_query = $blast{'query_name'};
  
  #Get info for exon feature
  $exon{'seqname'} = $blast{'query_name'};
  $exon{'start'} = 1;
  $exon{'end'} = $blast{'query_length'};
  $exon{'gene_id'} =  $blast{'hit_name'}; #This would be a very useful place to store cluster info
  $exon{'transcript_id'} =  $exon{'gene_id'} . '.1';
  if ( $blast{'hit_strand'} == 1 ) { $exon{'strand'} = '+'; }
  else { $exon{'strand'} = '-'; }
  print FlattenGTF(\%exon), "\n";
  
  #get info for cds feature
  $cds{'seqname'} = $exon{'seqname'};
  $cds{'start'} = $blast{'query_start'};
  $cds{'end'} = $blast{'query_end'};
  $cds{'score'} =  $blast{'score'}; 
  $cds{'gene_id'} =  $exon{'gene_id'};
  $cds{'strand'} = $exon{'strand'};
  $cds{'transcript_id'} =  $exon{'transcript_id'};
  print FlattenGTF(\%cds), "\n";
}   
