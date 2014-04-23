#!/usr/bin/perl -w

=head1 NAME

CompareSets.pl version 1, 24 June 2012

=head1 SYNOPSIS

CompareSets.pl file1 file2 [file3 file4] 

=head1 DESCRIPTION

Compares the values in the first column of three or four datasets and prints an R
command to generate a Venn Diagram of overlaps.

=head2 NOTES

use library(VennDiagram) in R to generate the plot

Options:

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-gmail-dot-comE<gt>

=cut
####################################################################################################

use strict;
use warnings;
use List::MoreUtils qw(uniq);
use File::Basename;
use List::Compare;

my @sets;
my @files; 
foreach(@ARGV) {
  #$_ =~ /(\w+)/;
  (my $name, my $path, my $in_ext) = fileparse($_, qr/\.[^.]*/);
  push(@files, $name);
  open(my $in, "<", $_);
  my @names;
  while(<$in>) {
    my @fields = split(/\s+/, lc($_));
    if ( @fields == 0 ) { next; }
    if ( $fields[0] =~ /\S/ ) {
      $fields[0] =~ s/^\s+//;     #strip out leading whitespace
      $fields[0] =~ s/\s+$//;     #strip out leading whitespace      
      push(@names, $fields[0]); }
  }
  @names = uniq(@names);
  push(@sets, \@names);
}

GetGroups(\@sets, \@files);
#GetVenn(\@sets, \@files);

sub GetGroups{
   my $set_ref = shift; 
   my @sets = @{$set_ref};
   my $file_ref = shift; 
   my @files = @{$file_ref};
   my @used = ();
   open(OUT, ">", 'groups.txt'); #This IO could be made a lot more flexible, but this will do for now
   my $lc= List::Compare->new(@sets);
   my @current = $lc->get_intersection;
   print join(", ", @files), ": ", scalar(@current), "\n";
   print OUT join(", ", @files), ": ", join(", ", @current), "\n";
   @used =(@used, @current);
   if ( @sets > 2 ) {
      print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
      for (my $i = 0; $i < @sets; $i ++) {
         my @spliced_set = @sets;
         my @spliced_files = @files;
         splice(@spliced_set, $i, 1);
         splice(@spliced_files, $i, 1);
         print join(", ", @spliced_files), ": ";
         print OUT join(", ", @spliced_files), ": ";
         $lc= List::Compare->new(@spliced_set);
         my @in_all = $lc->get_intersection;
         $lc= List::Compare->new(\@in_all, \@used);
         my @unique_to_set = $lc->get_unique;
         print OUT SplitArray(@unique_to_set), "\n";
         print scalar(@unique_to_set), "\n";
         @used =(@used, @unique_to_set);         
      }
   }
   if ( @sets > 3 ) {
      print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
      for (my $i = 0; $i < @sets - 1; $i ++) {
         for (my $j = $i + 1; $j < @sets; $j ++ ) {
            my @spliced_set = @sets;
            my @spliced_files = @files;
            splice(@spliced_set, $j, 1);
            splice(@spliced_set, $i, 1);
            splice(@spliced_files, $j, 1);
            splice(@spliced_files, $i, 1);
            $lc= List::Compare->new(@spliced_set);
            my @in_all = $lc->get_intersection;
            $lc= List::Compare->new(\@in_all, \@used);
            my @unique_to_set = $lc->get_unique;
            if ( @unique_to_set > 0 ) {
              print OUT join(", ", @spliced_files), ":\n";
              #print OUT join(", ", @unique_to_set), "\n";
              print OUT SplitArray(@unique_to_set), "\n";
            }
            print join(", ", @spliced_files), ": ";
            print scalar(@unique_to_set), "\n";
            @used =(@used, @unique_to_set);       
         }  
      }
   }
   if ( @sets > 4 ) {
      print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
      for (my $i = 0; $i < @sets - 2; $i ++) {
         for (my $j = $i + 1; $j < @sets; $j ++ ) {
            for ( my $k = $j + 1; $k < @sets; $k ++ ) {
               my @spliced_set = @sets;
               my @spliced_files = @files;
               splice(@spliced_set, $k, 1);
               splice(@spliced_set, $j, 1);
               splice(@spliced_set, $i, 1);
               splice(@spliced_files, $k, 1);
               splice(@spliced_files, $j, 1);
               splice(@spliced_files, $i, 1);
               $lc= List::Compare->new(@spliced_set);
               my @in_all = $lc->get_intersection;
               $lc= List::Compare->new(\@in_all, \@used);
               my @unique_to_set = $lc->get_unique;
               if ( @unique_to_set > 0 ) {
                 print OUT join(", ", @spliced_files), ":\n";
                 #print OUT join(", ", @unique_to_set), "\n";
                 print OUT SplitArray(@unique_to_set), "\n";
               }
               print join(", ", @spliced_files), ": ";
               print scalar(@unique_to_set), "\n";
               @used =(@used, @unique_to_set);    
            }   
         }  
      }
   }
   print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
   for (my $i = 0; $i < @sets; $i ++) {
      my @spliced_set = @{$sets[$i]}; 
      my $spliced_file = $files[$i];
      #print $spliced_file, ": ";
      #print OUT $spliced_file, ": ";
      $lc= List::Compare->new(\@spliced_set, \@used);
      my @unique_to_set = $lc->get_unique;
      if ( @unique_to_set > 0 ) {
         print "$spliced_file (", scalar(@unique_to_set), "):\n";
         #print OUT join(", ", @unique_to_set), "\n";
         #print OUT SplitArray(@unique_to_set), "\n";
         foreach (@unique_to_set) {
           print "$_\n";
         }
      }
      #print "$spliced_file: ";
      #print scalar(@unique_to_set), "\n";
      @used =(@used, @unique_to_set);       
   }
}

sub SplitArray {
  my @list = @_;
  my @split_array = ();
  while ( @list > 20 ) {
    print OUT @split_array, join(", ", @list[0..19]) . "\n";
    splice(@list, 0, 20);
  }
  if ( @list > 0 ) {
    print OUT @split_array, join(", ", @list) . "\n";
  }
  return @split_array;
}

sub GetVenn {
   my $set_ref = shift; 
   my @sets = @{$set_ref};
   my $file_ref = shift; 
   my @files = @{$file_ref};
   
   print "# drawing Venn Diagram for total of ", scalar(List::Compare->new(@sets)->get_union), " items (type library(VennDiagram) to view in R)\n";

   if ( @sets == 2 ) {
      print "ven.plot<-draw.pairwise.venn(\n";
   }
   elsif ( @sets == 3 ) {
      print "ven.plot<-draw.triple.venn(\n";
   }
   elsif ( @sets == 4 ) {
      print "venn.plot <- draw.quad.venn(\n";
   }
   elsif ( @sets == 5 ) {
      print "venn.plot <- draw.quintuple.venn(\n";
   }
   else { die "This script currently only produces pariwise to quintuple Venn plots\n"; }

   for (my $i = 0; $i < @sets; $i ++) { 
      print "  area", $i+1, " = ", scalar(@{$sets[$i]}), ",\n";
   }

   for (my $i = 0; $i < @sets - 1; $i ++) {
      for (my $j = $i + 1; $j < @sets; $j ++ ) {
         my $lc= List::Compare->new($sets[$i], $sets[$j]);
         if ( @sets == 2 ) {     
            print "  cross.area = ", scalar($lc->get_intersection), ",\n";
         }
         else {
            print "  n", $i+1, $j+1, " = ", scalar($lc->get_intersection), ",\n";
         }
      }
   }

   for (my $i = 0; $i < @sets - 2; $i ++) {
      for (my $j = $i + 1; $j < @sets - 1; $j ++ ) {
         for ( my $k = $j + 1; $k < @sets; $k ++ ) {
            my $lc= List::Compare->new($sets[$i], $sets[$j], $sets[$k]);
            print "  n", $i+1, $j+1, $k+1, " = ", scalar($lc->get_intersection), ",\n";
         }
      }
   }

   for (my $i = 0; $i < @sets - 3; $i ++) {
      for (my $j = $i + 1; $j < @sets - 2; $j ++ ) {
         for ( my $k = $j + 1; $k < @sets - 1; $k ++ ) {
            for ( my $l = $k + 1; $l < @sets; $l ++ ) {      
               my $lc= List::Compare->new($sets[$i], $sets[$j], $sets[$k], $sets[$l]);
               print "  n", $i+1, $j+1, $k+1, $l+1, " = ", scalar($lc->get_intersection), ",\n";
            }
         }
      }
   }

   for (my $i = 0; $i < @sets - 4; $i ++) {
      for (my $j = $i + 1; $j < @sets - 3; $j ++ ) {
         for ( my $k = $j + 1; $k < @sets - 2; $k ++ ) {
            for ( my $l = $k + 1; $l < @sets -1; $l ++ ) {
               for ( my $m = $l + 1; $m < @sets; $m ++ ) {    
                  my $lc= List::Compare->new($sets[$i], $sets[$j], $sets[$k], $sets[$l], $sets[$m]);
                  print "  n", $i+1, $j+1, $k+1, $l+1, $m+1, " = ", scalar($lc->get_intersection), ",\n";
               }
            }
         }
      }
   }

   print "  category = c('", join("','", @files), "'),\n";
   if ( @sets == 2) {
      print "  fill = c('blue', 'red'),
             cat.col = c('blue', 'red'),\n";
   }  
   elsif ( @sets == 3 ) {
      print "  fill = c('blue', 'red', 'green'),
         cat.col = c('blue', 'red', 'green'),\n";
   }
   elsif ( @sets == 4 ) {
      print "  fill = c('orange', 'red', 'green', 'blue'),
         cat.col = c('orange', 'red', 'green', 'blue'),\n";
   }
   else {
      print "  fill = c('orange', 'red', 'green', 'blue', 'purple'),
         cat.col = c('orange', 'red', 'green', 'blue', 'purple'),\n";
   }

   print "  lty = 'blank',
      cex = 2,
      cat.cex = 1.75,
      margin=0.05,
      fontfamily='sans',
      cat.fontfamily='sans'
      )\n";
}
#########################################################################################
open(my $out, ">", 'groups.txt');

