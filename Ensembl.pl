#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs(
{
    -HOST => 'mysql.ebi.ac.uk',
    -PORT => 4157,
    -USER => 'anonymous'
},
{
    -HOST => 'ensembldb.ensembl.org',
    -PORT => 5306,
    -USER => 'anonymous'
}
);

    
#my $slice_adaptor = $registry->get_adaptor( 'Selaginella moellendorffii', 'Core', 'Slice' );
#my $slice = $slice_adaptor->fetch_by_transcript_stable_id( 'EFJ31082', 1e2 );
#my $sequence = $slice->seq();
#print $sequence, "\n";

# get the MemberAdaptor
my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
    "plants","compara","Member");

# fetch a Member
my $id = 'EFJ35851';
my $member1 = $member_adaptor->fetch_by_source_stable_id(
    "ENSEMBLPEP",$id);

# print out some information about the Member
print $member1->chr_name, " ( ", $member1->chr_start, " - ", $member1->chr_end,
    " ): ", $member1->description, "\n";

my $gene = $member1->get_Gene;
my $member = $member_adaptor->fetch_by_source_stable_id(
    "ENSEMBLGENE",$gene->stable_id);

my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor("plants", "compara", "Homology");
my $homologies = $homology_adaptor->fetch_all_by_Member($member);

# That will return a reference to an array with all homologies (orthologues in
# other species and paralogues in the same one)
# Then for each homology, you can get all the Members implicated

foreach my $homology (@{$homologies}) {
  # You will find different kind of description
  # UBRH, MBRH, RHS, YoungParalogues
  # see ensembl-compara/docs/docs/schema_doc.html for more details
  print $homology->description," ", $homology->subtype, ": ";
  foreach my $member_attribute (@{$homology->get_all_Member_Attribute}) {
    my ($member, $attribute) = @{$member_attribute};
    unless ( $member->stable_id eq $id ) {
      print $member->taxon->genus, " ", $member->taxon->species, " ", $member->stable_id, " (", $attribute->perc_id, "%)","\n";
    }
  }
}


my $family_adaptor = Bio::EnsEMBL::Registry->get_adaptor("plants","compara","Family");
my $families = $family_adaptor->fetch_all_by_Member($member1);

foreach my $family (@{$families}) {
    print join(" ", map { $family->$_ }  qw(description description_score))."\n";

    foreach my $member_attribute (@{$family->get_all_Member_Attribute}) {
        my ($member, $attribute) = @{$member_attribute};
        print $member->stable_id," ",$member->taxon_id,"\n";
    }

    my $simple_align = $family->get_SimpleAlign();
    my $alignIO = Bio::AlignIO->newFh(
        -interleaved => 0,
        -fh          => \*STDOUT,
        -format      => "phylip",
        -idlength    => 20);

    print $alignIO $simple_align;

    $simple_align = $family->get_SimpleAlign("cdna");
    $alignIO = Bio::AlignIO->newFh(
        -interleaved => 0,
        -fh          => \*STDOUT,
        -format      => "phylip",
        -idlength    => 20);

    print $alignIO $simple_align;
}