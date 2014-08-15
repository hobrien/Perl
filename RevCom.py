#!/usr/local/bin/python

import sys
from Bio import SeqIO
if len(sys.argv) == 1:
  seqs = SeqIO.parse(sys.stdin, "fasta")
else:
  seqs = []
  for arg in sys.argv[1:]:
    seqs += SeqIO.parse(arg, "fasta")

revcom = []
for seq in seqs:
  print ">" + seq.id
  print seq.seq.reverse_complement()
