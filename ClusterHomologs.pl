#!/usr/bin/perl -w

=head1 NAME

ClusterHomologs.pl version 1, 7 May 2013

=head1 SYNOPSIS

cat BLASTRESULT | ClusterHomologs.pl > OUTFILE

=head1 DESCRIPTION

Clusters sequences by blast Evalue using single-linkage clustering


=head2 NOTES

Groups sequences that have significant blast hits. If a sequence has a significant blast
to any sequence in a group, then it is added to that group and if any sequences in
different groups have significant blast hits to each other the groups are merged.

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use warnings;
use strict;
#use Bio::Index::Fasta;
#use Bio::SeqIO;
use FileParser qw(ParseBlast FlattenBlast);
#$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";

my $usage = "cat BLASTFILE | TopHit.pl > OUTFILE";

my %clusters; #hash of arrays
my %cluster_id; #cluster membership for each sequence.
my %strand;
my $cluster_num = 0;

while (<>) {
  chomp;
  $_ =~ s/,\s*/\t/g;
  my %result = % {ParseBlast($_) };
  if ( $result{'query_name'} eq $result{'hit_name'} ) { next; }
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
foreach ( sort keys %cluster_id ) {
  print "$_\tCluster", "$cluster_id{$_}\t$strand{$_}\n";
}

#foreach ( sort {$a <=> $b} keys %clusters ) {
#  print "Cluster$_\t", scalar(@{$clusters{$_}}), "\t";
#  print join(",", @{$clusters{$_}}), "\n";
#}
    
