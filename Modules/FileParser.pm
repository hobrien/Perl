package FileParser;
require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK =qw(ParseBlast ParseGTF FlattenGTF FlattenBlast WriteFasta);
use strict;
use warnings;

sub WriteFasta {
  my $seq = shift;
  $seq =~ s/([^\n]{60})/$1\n/g;
  return $seq;
}

sub ParseGTF {
  my $line = shift;
  chomp($line);
  my @fields = split(/\t/, $line);
  my %tags;
  foreach (split(/;\s*/, pop(@fields) ) ) {
    $_ =~ s/"//g;
    my @attribute = split(/\s+/, $_);
    $tags{$attribute[0]} = $attribute[1];
  }
  $tags{'frame'} = pop(@fields);
  unless ( $tags{'frame'} =~ /^[012.]$/ ) { die "frame $tags{'frame'} not recognized. Must be 0, 1, 2 or .\n"; }
  $tags{'strand'} = pop(@fields);
  unless ( $tags{'strand'} =~ /^[+-.]$/ ) { die "strand $tags{'strand'} not recognized. Must be +, - or .\n"; }
  $tags{'score'} = pop(@fields);
  unless ( $tags{'score'} =~ /^[+-]?(?=\.?\d)\d*\.?\d*(?:e[+-]?\d+)?\z/i or $tags{'score'} eq '.') { die "score $tags{'score'} not recognized. Must be a number\n"; }
  $tags{'end'} = pop(@fields);
  unless ( $tags{'end'} =~ /^\d+$/ ) { die "end $tags{'end'} not recognized. Must be a positive integer"; }
  $tags{'start'} = pop(@fields);
  unless ( $tags{'start'} =~ /^\d+$/ ) { die "start $tags{'start'} not recognized. Must be a positive integer"; }
  $tags{'feature'} = pop(@fields);
  $tags{'source'} = pop(@fields);
  $tags{'seqname'} = pop(@fields);
  return \%tags;
}

sub FlattenGTF {
  my $ref = shift;
  my %tags = %{$ref};
  my @fields;
  unless (exists $tags{'seqname'} ) { die "no seqname tag\n"; }
  push(@fields, $tags{'seqname'} );
  delete  $tags{'seqname'};
  unless (exists $tags{'source'} ) { die "no source tag\n"; }
  push(@fields, $tags{'source'} ); 
  delete  $tags{'source'};
  unless (exists $tags{'feature'} ) { die "no feature tag\n"; }
  push(@fields, $tags{'feature'} ); 
  delete  $tags{'feature'};
  unless (exists $tags{'start'} ) { die "no start tag\n"; }
  push(@fields, $tags{'start'} ); 
  delete  $tags{'start'};
  unless (exists $tags{'end'} ) { die "no end tag\n"; }
  push(@fields, $tags{'end'} ); 
  delete  $tags{'end'};
  unless (exists $tags{'score'} ) { die "no score tag\n"; }
  push(@fields, $tags{'score'} ); 
  delete  $tags{'score'};
  unless (exists $tags{'strand'} ) { die "no strand tag\n"; }
  push(@fields, $tags{'strand'} ); 
  delete  $tags{'strand'};
  unless (exists $tags{'frame'} ) { die "no frame tag\n"; }
  push(@fields, $tags{'frame'} );
  delete  $tags{'frame'};
  my $attributes;
  foreach ( sort keys %tags ) { $attributes .= $_ . ' "' . $tags{$_} . '"; '; }
  return join("\t", (@fields, $attributes) );
}

sub ParseBlast {
  my $line = shift;
  chomp($line);
  my @fields = split(/\s+/, $line);
  #print scalar(@fields), "\n";
  my %tags;
  $tags{'query_name'} = shift(@fields);
  $tags{'query_length'} = shift(@fields);
  unless ( $tags{'query_length'} =~ /^\d+$/ ) { die "query_length $tags{'query_length'} not recognized. Must be a positive integer\n"; }
  $tags{'hit_name'} = shift(@fields);
  $tags{'hit_length'} = shift(@fields);
  unless ( $tags{'hit_length'} =~ /^\d+$/ ) { die "hit_length $tags{'hit_length'} not recognized. Must be a positive integer\n"; }
  $tags{'percent'} = shift(@fields);
  unless ( $tags{'percent'} =~ /^[+-]?(?=\.?\d)\d*\.?\d*(?:e[+-]?\d+)?\z/i ) { die "percent $tags{'percent'} not recognized. Must be a number\n"; }
  $tags{'aln_length'} = shift(@fields);
  unless ( $tags{'aln_length'} =~ /^\d+$/ ) { die "aln_length $tags{'aln_length'} not recognized. Must be a positive integer\n"; }
  $tags{'mismatch'} = shift(@fields);
  unless ( $tags{'mismatch'} =~ /^\d+$/ ) { die "mismatch count $tags{'mismatch'} not recognized. Must be a positive integer\n"; }
  $tags{'gap'} = shift(@fields);
  unless ( $tags{'gap'} =~ /^\d+$/ ) { die "gap count $tags{'gap'} not recognized. Must be a positive integer\n"; }
  my $start = shift(@fields);
  unless ( $start =~ /^\d+$/ ) { die "query start $start not recognized. Must be a positive integer\n"; }
  my $end = shift(@fields);
  unless ( $end =~ /^\d+$/ ) { die "query end $end not recognized. Must be a positive integer\n"; }
  if ( $start > $end ) {
    $tags{'query_strand'} = -1;
    $tags{'query_start'} = $end;
    $tags{'query_end'} = $start;
  }  
  else {
    $tags{'query_strand'} = 1;
    $tags{'query_start'} = $start;
    $tags{'query_end'} = $end;
  }
  $tags{'query_frame'} = shift(@fields);
  unless ( $tags{'query_frame'} =~ /^\-?[0123]$/ ) { die "query_frame $tags{'query_frame'} not recognized. Must be 1, 2, 3, -1, -2, -3 or 0\n"; }
  $start = shift(@fields);
  unless ( $start =~ /^\d+$/ ) { die "hit start $start not recognized. Must be a positive integer\n"; }
  $end = shift(@fields);
  unless ( $end =~ /^\d+$/ ) { die "hit end $end not recognized. Must be a positive integer\n"; }
  if ( $start > $end ) {
    $tags{'hit_strand'} = -1;
    $tags{'hit_start'} = $end;
    $tags{'hit_end'} = $start;
  }  
  else {
    $tags{'hit_strand'} = 1;
    $tags{'hit_start'} = $start;
    $tags{'hit_end'} = $end;
  }
  $tags{'hit_frame'} = shift(@fields);
  unless ( $tags{'hit_frame'} =~ /^\-?[0123]$/ ) { die "hit_frame $tags{'hit_frame'} not recognized. Must be 1, 2, 3, -1, -2, -3 or 0\n"; }
  $tags{'evalue'} = shift(@fields);
  unless ( $tags{'evalue'} =~ /^[+-]?(?=\.?\d)\d*\.?\d*(?:e[+-]?\d+)?\z/i ) { die "Evalue $tags{'evalue'} not recognized. Must be a number\n"; }
  $tags{'score'} = shift(@fields);
  unless ( $tags{'score'} =~ /^[+-]?(?=\.?\d)\d*\.?\d*(?:e[+-]?\d+)?\z/i ) { die "score $tags{'score'} not recognized. Must be a number\n"; }
  return \%tags;
}

sub FlattenBlast {
  my $ref = shift;
  my %tags = %{$ref};
  my @fields;
  unless (exists $tags{'query_name'} ) { die "no query_name tag\n"; }
  push(@fields, $tags{'query_name'} );
  delete  $tags{'query_name'};
  unless (exists $tags{'query_length'} ) { die "no query_length tag\n"; }
  push(@fields, $tags{'query_length'} ); 
  delete  $tags{'length_name'};
  unless (exists $tags{'hit_name'} ) { die "no hit_name tag\n"; }
  push(@fields, $tags{'hit_name'} ); 
  delete  $tags{'hit_name'};
  unless (exists $tags{'hit_length'} ) { die "no hit_length tag\n"; }
  push(@fields, $tags{'hit_length'} ); 
  delete  $tags{'hit_length'};
  unless (exists $tags{'percent'} ) { die "no percent tag\n"; }
  push(@fields, $tags{'percent'} ); 
  delete  $tags{'percent'};
  unless (exists $tags{'aln_length'} ) { die "no aln_length tag\n"; }
  push(@fields, $tags{'aln_length'} ); 
  delete  $tags{'aln_length'};
  unless (exists $tags{'mismatch'} ) { die "no mismatch tag\n"; }
  push(@fields, $tags{'mismatch'} ); 
  delete  $tags{'mismatch'};
  unless (exists $tags{'gap'} ) { die "no gap tag\n"; }
  push(@fields, $tags{'gap'} ); 
  delete  $tags{'gap'};
  unless (exists $tags{'query_start'} ) { die "no query_start tag\n"; }
  unless (exists $tags{'query_end'} ) { die "no query_end tag\n"; }
  unless (exists $tags{'query_strand'} ) { die "no query_strand tag\n"; }
  if ( $tags{'query_strand'} == 1 ) {
    push(@fields, ($tags{'query_start'},$tags{'query_end'}) ); 
  }
  elsif ( $tags{'query_strand'} == -1 ) {
    push(@fields, ($tags{'query_end'}, $tags{'query_start'}) ); 
  }
  else { 
    die "query_strand tag $tags{'query_strand'} not recognized. Should be 1 or -1\n"; 
  }
  delete  $tags{'query_start'};
  delete  $tags{'query_end'};
  delete  $tags{'query_strand'};
  unless (exists $tags{'query_frame'} ) { die "no query_frame tag\n"; }
  push(@fields, $tags{'query_frame'} );
  delete  $tags{'query_frame'};
  unless (exists $tags{'hit_start'} ) { die "no hit_start tag\n"; }
  unless (exists $tags{'hit_end'} ) { die "no hit_end tag\n"; }
  unless (exists $tags{'hit_strand'} ) { die "no hit_strand tag\n"; }
  if ( $tags{'hit_strand'} == 1 ) {
    push(@fields, ($tags{'hit_start'},$tags{'hit_end'}) ); 
  }
  elsif ( $tags{'hit_strand'} == -1 ) {
    push(@fields, ($tags{'hit_end'}, $tags{'hit_start'}) ); 
  }
  else { 
    die "hit_strand tag $tags{'hit_strand'} not recognized. Should be 1 or -1\n"; 
  }
  delete  $tags{'hit_start'};
  delete  $tags{'hit_end'};
  delete  $tags{'hit_strand'};
  unless (exists $tags{'hit_frame'} ) { die "no hit_frame tag\n"; }
  push(@fields, $tags{'hit_frame'} );
  delete  $tags{'hit_frame'};
  unless (exists $tags{'evalue'} ) { die "no evalue tag\n"; }
  push(@fields, $tags{'evalue'} );
  delete  $tags{'evalue'};
  unless (exists $tags{'score'} ) { die "no score tag\n"; }
  push(@fields, $tags{'score'} );
  delete  $tags{'score'};
  return join("\t", @fields);
}

1;