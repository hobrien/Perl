#!/usr/bin/perl -w

=head1 NAME

ConvertSeq.pl
2 Feb 2010, -modified 12 May 2010 to work on directories containing multiple files
modified 25 May 2011 to work with alignment files.

=head1 SYNOPSIS

ConvertSeqs.pl -i sequence file/folder -f outformat [options]

Options:
 -r replace existing outfile
 -x infile format (only necessary if it is not specified by the file extension)
 -n remove skip alignments with more or less than the specified # of sequences 
    (usually used to exclude alignments that do not include all taxa)
 -g remove gaps (skips alignment positions with a gap in any of the sequences
 -c concatinate alignments (this may give strange results if the alignments do not all have the same taxa)

=head1 DESCRIPTION
A simple script to convert sequence formats. 

If multiple sequences are present, each is written to the outfile

This seems to be working well on single sequences or folders of sequences.

=head2 NOTES

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################


use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::Utilities qw(:all);
use File::Basename;
use Getopt::Long;
use List::Util qw(min max);

my $infilename;
my $outfilename;
my $informat;
my $outformat;
my $help = 0;
my $replace = 0;
my $seq_num = 0;
my $remove_gaps = 0;
my $remove_strains = 0;
my $concat = 0;
my $remove_invariant = 0;
my $abb_name = 0;
my $remove_ambig = 0;
my $trunc_start = 0;
my $trunc_end = 0;
my @strains;

#Getopt::Long::Configure ("bundling");
GetOptions(
'infile:s' => \$infilename,
'format:s' => \$outformat,
'replace' => \$replace,
'xformat:s' => \$informat,
'number:s' => \$seq_num,
'outfile:s' => \$outfilename,
'gap' => \$remove_gaps,
'strain' => \$remove_strains,
'concatinate' => \$concat,
'variable' => \$remove_invariant,
'abbreviate' => \$abb_name,
'help|?' => \$help,
'ambiguous' =>\$remove_ambig,
'truncate:s' => \$trunc_start
);

my $usage = "type perldoc ConvertSeq.pl for help\n";
if( $help ) {
    die $usage;
}

if ( $trunc_start ) { 
  if ($remove_strains ) { die "Cannot select both truncate and remove seqs.\n"; }
  $trunc_end = shift; 
}
if ( $remove_strains) { 
  if ($trunc_start ) { die "Cannot select both truncate and remove seqs.\n"; }
  @strains = @ARGV; 
}
$infilename = shift unless $infilename;

$outformat = shift unless $outformat;
unless ($outformat) { die "output format not specified\n$usage";}
if ( $outformat eq "aln" ){ $outformat = "clustalw"; }
if ( $outformat eq "aln" ){ $outformat = "clustalw"; }
if ( $outformat eq "phy" ){ $outformat = "phylip"; }

####################################################################################################

(my $name, my $path, my $in_ext) = fileparse($infilename, qr/\.[^.]*/);

unless ($informat) { $informat = GetFormat($in_ext);}
if ( $informat eq "aln" ){ $informat = "clustalw"; }
my $out_ext = GetExtension($outformat);
unless ( $outfilename ) { $outfilename = $path . $name . $out_ext; }
if (-f $outfilename ) { 
  if ($replace) { print `rm $outfilename\n`; }
  else {die "outfile $outfilename already exists!\n"; }
}
if ( -d $infilename ) {
  opendir(DIR, $infilename) or die "can't opendir $infilename: $!";
  while (defined(my $filename = readdir(DIR))) {
    if ($filename =~ /^\./) { next; }
    $filename = "$infilename/$filename";
    ConvertSeqs($filename, $outfilename, $informat, $outformat);
  }
}

else {
    ConvertSeqs($infilename, $outfilename, $informat, $outformat);
}
exit;

####################################################################################################
sub ConvertSeqs {
  my $in_name = shift;
  my $out_name = shift;
  my $in_format = shift;
  my $out_format = shift;
  #output individual sequences (fasta, fastq, genbank)
  if ( $out_format eq "fasta" or $out_format eq "fastq" or $out_format eq "genbank") {
    my $outfile = Bio::SeqIO->new(-file => ">>$out_name" ,
                             '-format' => $out_format);
    #input individual sequences                         
    if ( $in_format eq "fasta" or $in_format eq "fastq" or $in_format eq "genbank") {
      my $infile = new Bio::SeqIO(-format => $in_format, 
                           -file   => $in_name);
      while ( my $seq = $infile->next_seq() ) {
        if ( $trunc_start ) { $seq = SeqTrunc($seq); }
        $outfile->write_seq($seq);
      }
    }
    elsif ( $in_format eq "snp_tbl" ) {
      my $aln = snp2aln($in_name);
      WriteSeqs($aln);
    }
    #input alignments (separate into individual sequences)
    else {
      my $infile = Bio::AlignIO->new(-format => $in_format, 
                           -file   => $in_name);
      while (my $aln = $infile->next_aln) {
        WriteSeqs($aln);
      } 
    } 
  }
#output is an alignment (input must be an alignment or snp table)
  else {
    if ( $in_format eq "snp_tbl" ) {
      my $aln = snp2aln($in_name);
      my $outfile;
      if ( $out_format =~ /phy/i and $out_format =~ /ext/i ) {  #extended phylip format
        open(OUT, ">>$out_name");
        WritePhylipExt($aln);
      }
      else { 
        $outfile = Bio::AlignIO->new(-file => ">>$out_name",
                                              '-format' => $out_format);
        $outfile->write_aln($aln);
      }
    }
    else {
      my $infile = Bio::AlignIO->new(-format => $in_format, 
                             -file   => $in_name);
      my $outfile;
      if ( $out_format =~ /phy/i and $out_format =~ /ext/i ) { open(OUT, ">>$out_name"); } #extended phylip format
      else {  $outfile = Bio::AlignIO->new(-file => ">>$out_name",
                                          '-format' => $out_format);}
      my $cat_aln;
      while ( my $aln = $infile->next_aln ) {
        if ($seq_num) {   #skip alignments without the specified number of sequences.
          my $num = scalar($aln->each_seq);
          unless ( $num == $seq_num ) { next; }
        }
#REMOVE GAPS (skipped if remove ambiguous also selected because the latter does both
        if ( $remove_gaps ) { unless ( $remove_ambig ) {$aln = $aln->remove_gaps; } }
#REMOVE AMBIGUOUS
        if ( $remove_ambig ) { $aln = RemoveAmbig($aln); }       
#REMOVE SELECTED STRAINS
        if ( $remove_strains ) {
          foreach ( @strains ) {
            my @seqs = $aln->each_seq_with_id($_);
            foreach (@seqs) { $aln->remove_seq($_); }
          }
        }
#REMOVE INVARIANT
        if ( $remove_invariant ) { $aln = $aln->remove_columns(['match']); }
#REMOVE ALIGNMENTS WITH NO SEQUENCE
        unless ( $aln->length ) { next; }
#ABBREVIATE NAMES
        if ( $abb_name ) { $aln = FixNames($aln); }
#TRUNCATE
        if ($trunc_start) { 
          foreach ( $aln->each_seq ) {
            $aln->remove_seq($_);
            my $seq = SeqTrunc($_);
            $aln->add_seq($seq);    
          }
        }
#CONCATINATE
        if ( $concat ) { 
          if ( $cat_aln ) { $cat_aln = cat($cat_aln, $aln); }
          else { $cat_aln = $aln; }
        }
        else {
         if ( $out_format =~ /phy/i and $out_format =~ /ext/i ) { WritePhylipExt($aln); }
         else { $outfile->write_aln($aln); }
        }
      }
      if ( $concat ) {
        if ( $out_format =~ /phy/i and $out_format =~ /ext/i ) { WritePhylipExt($cat_aln); }
        else { $outfile->write_aln($cat_aln); }
      }
    }
  }
}
######################################################################################
# IF ADDING NEW SEQUENCE FORMATS THEY MUST BE ADDED TO CONDITIONAL AT THE START OF THE ConvertSeqs SUBROUTINE
sub GetExtension {
  my $format = shift;
  if ( $format eq "fasta" ) { return ".fa";}
  elsif ( $format eq "xmfa" ) { return ".xmfa";}
  elsif ( $format eq "fastq" ) { return ".fq";}
  elsif ( $format =~ /phy/ ) { return ".phy";}
  elsif ( $format eq "clustalw" ) { return ".aln";}
  elsif ( $format eq "genbank" ) { return ".gbk";}
  elsif ( $format eq "nexus" ) { return ".nex";}
  else {die "file format not recognized!";}
}

sub GetFormat {
  my $ext = shift;
  my $format;
  if ( $ext =~ /f.*q/i ) { $format = "fastq";}
  elsif ( $ext =~ /^\.fa*/i ) { $format = "fasta";}
  elsif ( $ext =~ /xmfa/i ) { $format = "xmfa";}
  elsif ( $ext =~ /phy/i ) { $format = "phylip";}
  elsif ( $ext =~ /aln/i ) { $format = "clustalw";}
  elsif ( $ext =~ /gb/i ) { $format =  "genbank";}
  elsif ( $ext =~ /nxs/i ) { $format =  "nexus";}
  elsif ( $ext =~ /nex/i ) { $format =  "nexus";}
  elsif ( $ext =~ /maf/i ) { $format =  "maf";}
  elsif ( $ext =~ /txt/i ) { $format = "snp_tbl";}
  else {die "extension not recognized!";}
  return $format;
}

sub WritePhylipExt {
  my $aln = shift;
  my @seqs = $aln->each_seq;
  print OUT scalar(@seqs), " ", length($seqs[0]->seq), "\n";
  my @name_lengths;
  foreach (@seqs) { push(@name_lengths, length($_->display_id)); }
  my $name_size = max(@name_lengths) + 1;
  if ( $name_size < 10 ) { $name_size = 10; }
  foreach (@seqs) {
    printf OUT "%-*s",  $name_size, $_->display_id;
    print OUT $_->seq, "\n";
  }
}

sub FixNames {
  my $aln = shift;
  foreach ( $aln->each_seq ) {
    $aln->remove_seq($_);
    $_->display_id =~ /\/.*\/([^.]+)/;
    $_->display_id($1);
    $aln->add_seq($_);
  }
  return $aln;
}

sub RemoveAmbig {
  my $aln = shift;
  foreach ( $aln->each_seq ) {
    $aln->remove_seq($_);
    my $seq = $_->seq;
    $seq =~ s/n/-/gi;
    $_->seq($seq);
    $aln->add_seq($_);
  }
  $aln = $aln->remove_gaps;
  return $aln;
}

sub SeqTrunc {
  my $seq = shift;
  unless ( $trunc_end ) { $trunc_end = $seq->length; }
  return $seq->trunc($trunc_start, $trunc_end);
}

sub WriteSeqs {
  my $aln = shift;
  my $outfile = Bio::SeqIO->new('-file' => ">$outfilename",
         '-format' => $outformat) or die "could not open seq file $outfilename\n\n$usage";
  foreach ( $aln->each_seq) {
    my $seq = $_;
    if ( $trunc_start ) { $seq = SeqTrunc($seq); }
    my $string = $seq->seq;
    $string =~ s/-//g;
    $seq->seq($string);
    $outfile->write_seq($seq);
  }
}

sub snp2aln {
  my $file = shift;
  open (IN, $file);
  my @taxa;
  my @seqs;
  while (<IN>) {
    chomp($_);
    $_ =~ s/\s+$//;
    my @chars = split(/\t/, $_);
    shift(@chars);
    shift(@chars);
    unless (@taxa) { @taxa = @chars; next; }
    for ( my $x = 0; $x < @chars; $x ++ ) { $chars[$x] =~ s/.*\(([ACGT])\)/$1/; }
    unless (@seqs) { @seqs = @chars; next; }
    for ( my $x = 0; $x < scalar(@chars); $x ++ ) {
      $seqs[$x] .= $chars[$x]
    }
  }
  my $aln = new Bio::SimpleAlign;
  for ( my $x = 0; $x < scalar(@seqs); $x ++ ) {
    my $seq = new Bio::LocatableSeq(-seq => $seqs[$x],
                                    -id => $taxa[$x],
                                    -start => 1,
                                    -end => length($seqs[$x]));
    $aln->add_seq($seq);
  }
  return $aln;
}