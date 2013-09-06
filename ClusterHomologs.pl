#!/usr/bin/perl -w

=head1 NAME

ClusterHomologs.pl version 1, 7 May 2013

=head1 SYNOPSIS

cat BLASTRESULT | ClusterHomologs.pl evalue > OUTFILE

=head1 DESCRIPTION

Clusters sequences by blast Evalue using single-linkage clustering

=head2 NOTES

Groups sequences that have significant blast hits. If a sequence has a significant blast
to any sequence in a group, then it is added to that group and if any sequences in
different groups have significant blast hits to each other the groups are merged.

Generates 2 outfiles: 
  cluster_summary.txt: info about cluster size and membership for each cluster
  cluster_id: cluster number and orientation for each sequence

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use warnings;
use strict;
use FileParser qw(ParseBlast);
use List::Util qw(min max);

my $usage = "cat BLASTFILE | ClusterHomologs.pl evalue";

my %clusters; #hash of arrays
my %cluster_id; #cluster membership for each sequence.
my %starts; #minimum start position of sequence in blast results
my %ends; #maximum end position of sequence in blast results
my %strand;
my $cluster_num = 0;
my $evalue = shift;
unless ( $evalue ) { $evalue = 10; }

while (<>) {
  chomp;
  $_ =~ s/,\s*/\t/g;
  my %result = % {ParseBlast($_) };
  if ( $result{'query_name'} eq $result{'hit_name'} or $result{'evalue'} > $evalue ) { next; }
  unless ($starts{$result{'query_name'}} and $starts{$result{'query_name'}} < $result{'query_start'} ) { 
    $starts{$result{'query_name'}} = $result{'query_start'};
  }
  unless ($ends{$result{'query_name'}} and $ends{$result{'query_name'}} > $result{'query_end'} ) { 
    $ends{$result{'query_name'}} = $result{'query_end'};
  }
  unless ($starts{$result{'hit_name'}} and $starts{$result{'hit_name'}} < $result{'hit_start'} ) { 
    $starts{$result{'hit_name'}} = $result{'hit_start'};
  }
  unless ($ends{$result{'hit_name'}} and $ends{$result{'hit_name'}} > $result{'hit_end'} ) { 
    $ends{$result{'hit_name'}} = $result{'hit_end'};
  }
  if ( $cluster_id{$result{'query_name'}} ) {
    if ( $cluster_id{$result{'hit_name'}} ) {
      unless ( $cluster_id{$result{'query_name'}} == $cluster_id{$result{'hit_name'}}) { #sequences are part of different clusters (merge)
        my $old_cluster = $cluster_id{$result{'hit_name'}};
        foreach ( @{$clusters{$cluster_id{$result{'hit_name'}}}} ) { 
          $cluster_id{$_} = $cluster_id{$result{'query_name'}};
          push(@{$clusters{$cluster_id{$result{'query_name'}}}}, $_);
          #modify strand information
          if ( $result{'hit_strand'} == -1 ) { $strand{$_} = 0 - $strand{$_}; }   
        }
        delete $clusters{$old_cluster};
      }
    }  # if hit and query sequences are already in the same cluster, nothing more to be done.
    else { #add hit sequence to query cluster
      $cluster_id{$result{'hit_name'}} = $cluster_id{$result{'query_name'}};
      push(@{$clusters{$cluster_id{$result{'query_name'}}}}, $result{'hit_name'});
    #record strand information
    if ( $result{'hit_strand'} == 1 ) { $strand{$result{'hit_name'}} = $strand{$result{'query_name'}}; }
    else { $strand{$result{'hit_name'}} = 0 - $strand{$result{'query_name'}}; }   
    }
  }
  elsif ( $cluster_id{$result{'hit_name'}} ) {  #add query sequence to hit cluster
    $cluster_id{$result{'query_name'}} = $cluster_id{$result{'hit_name'}};
    push(@{$clusters{$cluster_id{$result{'hit_name'}}}}, $result{'query_name'});
    #record strand information
    if ( $result{'hit_strand'} == 1 ) { $strand{$result{'query_name'}} = $strand{$result{'hit_name'}}; }
    else { $strand{$result{'query_name'}} = 0 - $strand{$result{'hit_name'}}; }   
  }
  else { #hit and query constitute a new cluster
    $cluster_num ++;
    $cluster_id{$result{'query_name'}} = $cluster_num; 
    $cluster_id{$result{'hit_name'}} = $cluster_num;
    $clusters{$cluster_num} = [$result{'query_name'}, $result{'hit_name'} ];
    #record strand information
    $strand{$result{'query_name'}} = 1;
    if ( $result{'hit_strand'} == 1 ) { $strand{$result{'hit_name'}} = 1; }
    else { $strand{$result{'hit_name'}} = -1; }
  }
}

open (my $sum, ">", 'cluster_summary.txt');
foreach my $cluster ( sort {$a <=> $b} keys %clusters ) {
  print $sum "Cluster$cluster\t", scalar(@{$clusters{$cluster}}), "\t";
  print $sum join(",", @{$clusters{$cluster}}), "\n";
  my $ref_strand;
  foreach (@{$clusters{$cluster}}) {
    if ( $_ =~ /EFJ/i or $_ =~ /ADH/i or $_ =~ /jgi/i or $_ =~ /SELMO/i ) { #reference transcripts
      if ( $ref_strand and $ref_strand != $strand{$_} ) {
        warn "Cluster $cluster contains reference sequences on opposite strands\n";
        $ref_strand = 1;
      }
      else { $ref_strand = $strand{$_}; }
    }
  }
  if ( $ref_strand and $ref_strand == -1 ) {
    foreach (@{$clusters{$cluster}}) { $strand{$_} = 0 - $strand{$_}; }
  }
}
close ($sum);

open(my $seq_info, ">", 'cluster_ids.txt');
foreach ( sort keys %cluster_id ) {
  print $seq_info "$_\tCluster", "$cluster_id{$_}\t$strand{$_}\t$starts{$_}\t$ends{$_}\n";
  #print  "$_\tCluster", "$cluster_id{$_}\t$strand{$_}\t$starts{$_}\t$ends{$_}\n";
}

    
