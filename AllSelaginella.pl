#!/usr/bin/perl -w

=head1 NAME

AllSelaginella.pl version 1, 24 April 2013

=head1 SYNOPSIS

AllSelaginella.pl 'COMMAND XXX#.fa'

=head1 DESCRIPTION

Executes that command between quotes on all Selaginella datasets. XXX is replaced by the
species name (KRAUS/MOEL/UNC/WILD) and # is replaced by the sample number (1/2/3/4).


=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut


####################################################################################################

use warnings;
use strict;

my $command = shift;

foreach my $species ( qw(KRAUS MOEL WILD UNC) ) {
  foreach my $num ( (1,2,3,4) ) {
    if ( $species eq 'MOEL' and $num == 4 ) { next; }
    my $com = $command;
    $com =~ s/XXX/$species/g;
    $com =~ s/#/$num/g;
    print `$com`;
  }
}