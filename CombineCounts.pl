#!/usr/bin/perl -w

=head1 NAME

CombineCounts.pl version 1, 23 June 2012

=head1 SYNOPSIS

CombineCounts.pl > outfile 

=head1 DESCRIPTION

Combines count output from MOEL_DESeq.R, KRAUS_DESeq.R, UNC_DESeq.R and WILD_DESeq.R,
providing counts for each gene in each species and indicating values that are significant 
in at least one comparison. 

=head2 NOTES

A better way to do things would be to combine the raw counts and run DESeq on these data
so that comparisons can be made across species in absolute counts.

Options:

None. All file names are hard coded.

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use warnings;
use strict;
use autodie qw(open close);
use List::Util qw(sum);
use File::HomeDir;
use FileParser qw(ParseGTF);

my $HOME = File::HomeDir->my_home;
my ( $count_ref, $multi_ref ) = GetCounts('MOEL');
my %MOEL_counts = % { $count_ref };
my %MOEL_multi = % { $multi_ref };
( $count_ref, $multi_ref ) = GetCounts('WILD');
my %WILD_counts = % { $count_ref };
my %WILD_multi = % { $multi_ref };
( $count_ref, $multi_ref ) = GetCounts('KRAUS');
my %KRAUS_counts = % { $count_ref };
my %KRAUS_multi = % { $multi_ref };
( $count_ref, $multi_ref ) = GetCounts('UNC');
my %UNC_counts = % { $count_ref };
my %UNC_multi = % { $multi_ref };

print "Gene,MOEL1,MOEL2,MOEL3,MOEL_SUM,UNC1,UNC2,UNC3,UNC4,UNC_sum,WILD1,WILD2,WILD3,WILD4,WILD_SUM,KRAUS1,KRAUS2,KRAUS3,KRAUS4,KRAUS_SUM\n";
foreach ( sort keys %MOEL_counts ) { 
  if ( $MOEL_multi{$_} ) { warn "$_: multiple gene_ids in MOEL\n"; }
  print "$_,$MOEL_counts{$_}";
  delete $MOEL_counts{$_};
  if ( $UNC_multi{$_} ) { warn "$_: multiple gene_ids in UNC\n"; }
  if ( $UNC_counts{$_} ) { 
    print ",$UNC_counts{$_}"; 
  }
  else { print ",NA,NA,NA,NA,NA"; }
  delete $UNC_counts{$_};
  if ( $WILD_multi{$_} ) { warn "$_: multiple gene_ids in WILD\n"; }
  if ( $WILD_counts{$_} ) { 
    print ",$WILD_counts{$_}"; 
  }
  else { print ",NA,NA,NA,NA,NA"; }
  delete $WILD_counts{$_};
  if ( $KRAUS_multi{$_} ) { warn "$_: multiple gene_ids in KRAUS\n"; }
  if ( $KRAUS_counts{$_} ) { 
    print ",$KRAUS_counts{$_}\n"; 
  }
  else { print ",NA,NA,NA,NA,NA\n"; }
  delete $KRAUS_counts{$_};
}

foreach ( sort keys %UNC_counts ) { 
  if ( $UNC_multi{$_} ) { warn "$_: multiple gene_ids in UNC\n"; }
  print "$_,NA,NA,NA,NA,$UNC_counts{$_}";
  delete $UNC_counts{$_};
  if ( $WILD_multi{$_} ) { warn "$_: multiple gene_ids in WILD\n"; }
  if ( $WILD_counts{$_} ) { 
      print ",$WILD_counts{$_}"; 
  }
  else { print ",NA,NA,NA,NA,NA"; }
  delete $WILD_counts{$_};
  if ( $KRAUS_multi{$_} ) { warn "$_: multiple gene_ids in KRAUS\n"; }
  if ( $KRAUS_counts{$_} ) { 
    print ",$KRAUS_counts{$_}\n"; 
  }
  else { print ",NA,NA,NA,NA,NA\n"; }
  delete $KRAUS_counts{$_};
}

foreach ( sort keys %WILD_counts ) { 
  if ( $WILD_multi{$_} ) { warn "$_: multiple gene_ids in WILD\n"; }
  print "$_,NA,NA,NA,NA,NA,NA,NA,NA,NA,$WILD_counts{$_}";
  delete $WILD_counts{$_};
  if ( $KRAUS_multi{$_} ) { warn "$_: multiple gene_ids in KRAUS\n"; }
  if ( $KRAUS_counts{$_} ) { 
    print ",$KRAUS_counts{$_}\n"; 
  }
  else { print ",NA,NA,NA,NA,NA\n"; }
  delete $KRAUS_counts{$_};
}

foreach ( sort keys %KRAUS_counts ) { 
  if ( $KRAUS_multi{$_} ) { warn "$_: multiple gene_ids in KRAUS\n"; }
  print "$_,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,$KRAUS_counts{$_}\n";
  delete $KRAUS_counts{$_};
}


sub GetCounts {
  my $name = shift;
  my $gtffilename = "$HOME/Bioinformatics/Selaginella/GTF/" . $name . "_clc_scaf.gtf";
  open (my $gtf, "<", $gtffilename);
  my %gene_names;
  while (<$gtf>){
    my %tags = % { ParseGTF($_) };
    $gene_names{$tags{'gene_id'}} = $tags{'gene_name'};  
  }
  
  my %sig1;
  my %sig2;
  my %sig3;
  my %sig4;
  my @sig = ( \%sig1, \%sig2, \%sig3, \%sig4 );
  
  my $pvalfile = "$HOME/Documents/Selaginella/" . $name . "_PVAL.csv";
  open(my $pval, '<', $pvalfile);
  while(<$pval>){
    if ( $. == 1 ) { next; }
    chomp;
    $_ =~ s/"//g;
    my @fields = split(/,/, $_);
    if ( $fields[-2] eq 'NA' or $fields[-2] > 0.1 ) { next; }
    for ( my $i = 1; $i < 5; $i ++ ) {
      if ( $fields[-1] =~ /$i/ ) { $sig[$i-1]{$fields[1]} = 1; }  
    }
  }

  my %counts;
  my %multiple;
  my $countfilename = "$HOME/Documents/Selaginella/" . $name . "_ALL.csv";
  open(my $countfile, '<', $countfilename);
  while(<$countfile>) {
    if ( $. == 1 ) { next; }
    chomp;
    $_ =~ s/"//g;
    my @fields = split(/,/, $_);
    my $gene_id = shift(@fields);
    my $total = sum(@fields);
    push(@fields, $total);
    for ( my $i = 0; $i < @fields; $i ++ ) {
      if ( $fields[$i] < 10 ) {$fields[$i] = sprintf "%.3f", $fields[$i]; }
      elsif ( $fields[$i] < 100 ) {$fields[$i] = sprintf "%.2f", $fields[$i]; }
      elsif ( $fields[$i] < 1000 ) {$fields[$i] = sprintf "%.1f", $fields[$i]; }
      else { $fields[$i] = int($fields[$i]); }
      if ( $i == @fields - 1 ) { next; }
      if ( $sig[$i]{$gene_id} ) { $fields[$i] .= '*'; }
    }
    if ( $gene_names{$gene_id} ) { 
      if ( $counts{$gene_names{$gene_id}} ) {  #multiple genes match the same ref. Need to pick the one with the highest total count unless the one with the lower count is significantly DE
        $multiple{$gene_names{$gene_id}} = 1;
        $counts{$gene_names{$gene_id}} =~/([^,])$/;
        if ( $1 < $fields[-1] ) {
          unless ( $counts{$gene_names{$gene_id}} =~ /\*/ and join(",", @fields) !~ /\*/ ) { 
            $counts{$gene_names{$gene_id}} = join(",", @fields);
          }
        }
        elsif ( $counts{$gene_names{$gene_id}} !~ /\*/ and join(",", @fields) =~ /\*/ ) {
           $counts{$gene_names{$gene_id}} = join(",", @fields);
        } 
      }
      else { $counts{$gene_names{$gene_id}} = join(",", @fields); }
    }
  }
  return (\%counts, \%multiple);
}

