#!/opt/local/bin/python

import sys
from os import system, path

def main(argv):
  fasta = argv[0]
  phylip = fasta.replace('.fa', '.phy')
  phyml_tree = phylip + '_phyml_tree.txt'
  phyml_stats = phylip + '_phyml_stats.txt'
  tree = fasta.replace('.fa', '.nwk')
  if not path.exists(phylip):
    system("mafft --quiet %s |ConvertAln.pl -x fasta -f phyext -o %s -r" % (fasta, phylip))
  if not path.exists(tree):
    system("phyml --quiet --no_memory_check -i %s" % phylip)
    system("mv %s %s" % (phyml_tree, tree))
    system("rm %s" % phyml_stats)

if __name__ == "__main__":
   main(sys.argv[1:])
