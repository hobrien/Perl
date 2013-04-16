

my $file1name = shift;
my $index1 = shift;
$index1 --;
my $file2name = shift;
my $index2 = shift;
$index2 --;
open($file1, $file1name)

my %file1_genes;
while (<$file1>) {
  if ( $. ==0 ) { next; }
  chomp;
  my @fields = split(/\s+/, $_);
  $file1_genes{$fields[$index]} = 1;
}

