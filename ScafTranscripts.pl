#!/usr/bin/perl -w

=head1 NAME

ScafTranscripts.pl version 1, 2 April 2013

=head1 SYNOPSIS

ScafTranscripts.pl -s contig_seqs.fa -r reference_seqs.fa -b blast_results.bl [-o out_file] [-l logfile.txt]

=head1 DESCRIPTION

=head2 NOTES


=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-utoronto-dot-caE<gt>

=cut
####################################################################################################
use warnings;
use strict;
use Bio::Index::Fasta;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;
use FileParser qw(ParseBlast FlattenBlast);
use autodie qw(open close);
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";

my $blastfilename;
my $contigfilename;
my $reffilename;
my $logfilename;
my $outfilename;
GetOptions( 'seqfilename=s' => \$contigfilename,
	        'blastfilename=s' => \$blastfilename,
	        'outfilename=s' => \$outfilename,
	        'logfilename=s' => \$logfilename
	      )
  or die "failure to communicate\n";

#Prepare contig sequence index
my $contig_inx_name = $contigfilename . ".inx";
my $contig_inx = Bio::Index::Fasta->new(-filename => $contig_inx_name, -write_flag => 1);
#unless (-f $contig_inx_name ) {
  $contig_inx->make_index($contigfilename);
#}


#Open infle, outfile and log file
(my $name, my $path, my $ext) = fileparse($contigfilename, qr/\.[^.]*/);
unless ( $outfilename ) { $outfilename = $path . $name . '_scaf.fa'; }
unless ( $logfilename ) { $logfilename = $name . '_scaf_log.csv'; }
my $outfile = Bio::SeqIO->new('-file' => ">$outfilename",
         '-format' => 'fasta') or die "could not open seq file $outfilename\n";
open (my $log, ">", $logfilename);
#print $log "Query_ID1,Query_length1,Hit_ID1,Hit_length1,Percent_ID1,Aligned_length1,Mismatch_count1,Gap_count1,Query_start1,Query_end1,Hit_start1,Hit_end1,Evalue1,Score1,Num_hits1,Multi_hits1,Orientation1,Gap_size,Extra_seq,Query_ID2,Query_length2,Hit_ID2,Hit_length2,Percent_ID2,Aligned_length2,Mismatch_count2,Gap_count2,Query_start2,Query_end2,Hit_start2,Hit_end2,Evalue2,Score2,Num_hits2,Multi_hits2,Orientation2,Scaffold\n";
open (my $infile, "<", $blastfilename);

# make a hash of arrays containing the info of interest, with hit names as keys
my %hits;
my $previous = 'init';
while (<$infile>) {
  #print $., "\n";
  #next unless $. > 1;   #skip the first line
  my %result = %{ ParseBlast($_) };
  if ( $result{'query_name'} eq $previous ) { next; } #skip all but first line for each query_id
  $previous = $result{'query_name'};
  if ( $result{'percent'} < 70 ) { next; }
  my $query_seq = $contig_inx->fetch($result{'query_name'})  or die "could not find sequence $result{'query_name'}\n";
  my $skip = 0;
  if ( $result{'hit_start'} > 15  and $result{'query_start'} > 45 ) { $skip =1;}   #skip contigs with >45 bp of non-homologous overlap with hit
  if ( $result{'hit_length'} - $result{'hit_end'} > 15 and $result{'query_length'} - $result{'query_end'} > 45 ) {$skip = 1; }
  if ( $skip ) {
    if ( $result{'hit_strand'} == -1 ) { $query_seq = $query_seq->revcom; }
    $outfile->write_seq($query_seq);
    print $log FlattenBlast(\%result), "\n";
  }
  elsif ( $hits{$result{'hit_name'}} ) {
    my @contigs = @{$hits{$result{'hit_name'}}};
    my $index = 0;
    foreach ( @contigs ) {
      my %contig_info = %{ ParseBlast($_) };
      if ( $contig_info{'hit_start'} < $result{'hit_start'} ) { $index ++; }
      else { last; }
    }
    splice(@contigs, $index, 0, FlattenBlast(\%result));
    $hits{$result{'hit_name'}} = \@contigs;
  }
  else { 
    my @new_contig = ( FlattenBlast(\%result) );
    $hits{$result{'hit_name'}} = \@new_contig;
  }
}
close($infile);

my $scaf_num = 1;
foreach (keys %hits ) {
  my $contig_num = 1;
  my @contigs = @{$hits{$_}};
  my $scaf_string;
  my $scaf_name = $name . '_scaffold_' . $scaf_num;
  if ( scalar(@contigs) > 1 ) { 
    for ( my $x = 1; $x < scalar(@contigs); $x ++ ) {
      my %contig1 = %{ ParseBlast($contigs[$x-1]) };
      my $seq1 = $contig_inx->fetch($contig1{'query_name'})  or die "could not find sequence $contig1{'query_name'}\n";
      if ( $contig1{'query_strand'} == -1 ) { $seq1 = $seq1->revcom; }
      unless ( $scaf_string ) { $scaf_string = $seq1->seq; }
      my %contig2 = % { ParseBlast($contigs[$x]) };
      my $seq2 = $contig_inx->fetch($contig2{'query_name'})  or die "could not find sequence $contig2{'query_name'}\n";
      if ( $contig2{'query_strand'} == -1 ) { $seq2 = $seq2->revcom; }
      my $seq2string = $seq2->seq;
      if ( $contig2{'hit_end'} <= $contig1{'hit_end'} ) {
        print $log FlattenBlast(\%contig2), "\tredundant to $contig1{'query_name'}\n";
        splice(@contigs,$x, 1);
        $x --;
      }
      else {
        $contig_num ++;
        my $gap = ($contig2{'hit_start'} - $contig1{'hit_start'} -1) * 3;
        my $overhang1 = $contig1{'query_length'} - $contig1{'query_end'};
        my $overhang2 = $contig2{'query_start'};
        if ( $gap < 0 ) {  #contigs overlap; trim first to end of hit and second to
          $scaf_string = substr($scaf_string, 0, length($scaf_string) - $contig1{'query_length'} + $contig1{'query_end'}); #trim non-homologous sequence form first seq
          $scaf_string .= substr($seq2string, $contig2{'query_start'} - $gap - 1); #trim overlap + seq before hit from second seq
        }
        else {
          my $num_Ns = $gap - $overhang1 - $overhang2;
          my $trim1 = 0;
          my $trim2 = 0;
          while ( $num_Ns < 0 ) {
            if ($trim1 == $overhang1 and $trim2 == $overhang2 ) {
              print "This shouldn't be possible, but I can't think of anywhere else that there could be a infinite loop\n";
              die;
            }
            if ( $overhang1 > $trim1 ) { $trim1 ++; $num_Ns ++; }
            if ( $overhang2 > $trim2 ) { $trim2 ++; $num_Ns ++; }
          }
          $scaf_string = substr($scaf_string, 0, length($scaf_string) - $trim1) . 'N' x $num_Ns . substr($seq2string, $trim2);
        }
        my $extra = $overhang1 + $overhang2;
        if ( $contig_num == 2 ) {
          print $log FlattenBlast(\%contig1), "\t$scaf_name\tstart\n";
        }
        print $log FlattenBlast(\%contig2), "\t$scaf_name\tGap: $gap\tNon-homologous bp: $extra\n";
      }
    }
  }
  if ( $contig_num > 1 ) {
    $scaf_num ++;
    my $scaf_seq = Bio::Seq->new(-seq => $scaf_string,
                                 -alphabet => 'dna',
                                 -id  => $scaf_name);
    $outfile->write_seq($scaf_seq);  
  }
  else { 
    my %contig = %{ ParseBlast($contigs[0]) };
    my $seq = $contig_inx->fetch($contig{'query_name'})  or die "could not find sequence $contig{'query_name'}\n";
    if ( $contig{'query_strand'} == -1 ) { $seq = $seq->revcom; }
    $outfile->write_seq($seq);   
    print $log FlattenBlast(\%contig), "\n";  
  }
}
close($log);
$outfile->close;
exit;
