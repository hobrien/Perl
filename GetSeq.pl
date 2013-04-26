#!/usr/bin/perl -w

=head1 NAME

GetSeq.pl version 2, 23 April 2013

=head1 SYNOPSIS

ContigSort.pl [-n acc] file_name sequence_name

=head1 DESCRIPTION

Indexes fasta file and retrieves specified sequences (either individual name or a stream from STDIN)

=head2 NOTES

options:
 -n acc will cause the indexing to parse NCBI fasta headers to extract the accession number
 
=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Index::Fasta;
$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";

use Getopt::Long;

#Read in command line arguments
my $format;

Getopt::Long::Configure ("bundling");
GetOptions(
'n|name:s' => \$format,

);

my $infilename = shift;
my $inx_name = $infilename . ".inx";
my $inx = Bio::Index::Fasta->new(-filename => $inx_name, -write_flag => 1);

if ( $name and $name =~ /^acc/ ) {
  $inx->id_parser(\&get_acc);
}

$inx->make_index($infilename);

my $name = shift;

my $seqout = Bio::SeqIO->new(
                             -file => "| cat",
                             -format => 'Fasta',
                             );

if ( $name ) {
  my $seq = $inx->fetch($name) or die "could not find sequence $name\n";
  $seqout->write_seq($seq);
}
else {
  while (<>){
    chomp;
    my $seq = $inx->fetch($_) or die "could not find sequence $_\n";

  $seqout->write_seq($seq);
  }
}


sub get_acc {
   my $header = shift;
   $header =~ /^>gi\|\d+\|\w+\|(\w+)\.\d\|/;
   $1;
}