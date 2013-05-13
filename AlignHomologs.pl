#!/usr/bin/perl -w

=head1 NAME

AlignHomologs.pl version 1, 7 May 2013

=head1 SYNOPSIS

cat Cluster_Names | AlignHomologs.pl

=head1 DESCRIPTION

Aligns sequences in each sequence file made by GroupSequences.pl, trims ambiguous regions and makes trees


=head2 NOTES

Currently, this only alignes sequences from clusters that are represented by all species

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################


while (<>){
  if ( $_ =~ /(EFJ)|(ADH)|(jgi)/ and $_ =~ /UNC/ and $_ =~ /MOEL/ and $_ =~ /WILD/ and $_ =~ /KRAUS/ ) {
    $_ =~ /(\w+)/;
    print "Aligning $1\n";
    $_ =~ /(Cluster\d+)/;
    print `mafft --quiet Clusters/$1.fa >Alignments/$1.aln.fa`;
    print `trimal -in Alignments/$1.aln.fa -out Trimmed/$1.trim.fa -gappyout`;
    print `ConvertSeq.pl -i Trimmed/$1.trim.fa -f phyext -o Phylip/$1.phy`;
    print `phyml -i Phylip/$1.phy`;
    print `mv Phylip/$1.phy_phyml_tree.txt Trees/$1.nwk`;
  }
}