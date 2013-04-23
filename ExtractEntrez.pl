#!/usr/bin/perl -w

=head1 NAME

ExtractEntrez.pl -v 1 22 April 2013

=head1 SYNOPSIS

cat seq_file | ExtractEntrez.pl > out_file 

=head1 DESCRIPTION
Searches for genbank accession numbers in a file (sequences or list) and prints
an entrez query that can be used to exclude these sequences from the output of a blast
search. Useful for searching for new sequences without returning ones that are already
in your dataset

=head2 NOTES

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################


my @accessions;
while (<>){
  chomp;
  if ( $_ =~ /gi\|\d+\|\w+\|(\w+)/ ) { push(@accessions, $1); }
}

my $prefix = 'init';
my $first = 0;
my $previous = 0;

foreach (sort @accessions) {
  $_ =~ /(\D+)(\d+)/;
  if ( $1 eq $prefix and $2 == $previous + 1 ) { $previous ++; next; }
  unless ( $prefix eq 'init' ) {
    print ' NOT ', $prefix, $first;
    unless ( $first == $previous ) { print ':', $prefix, $previous; }
    print '[ACCN]';
  }
  $prefix = $1;
  $first = $2;
  $previous = $2;
}

print ' NOT ', $prefix, $first;
unless ( $first == $previous ) { print ':', $previous; }
print "[ACCN]\n";
