#!/opt/local/bin/python

import sys, getopt
from os import listdir, path, system
from Bio import SeqIO
from StringIO import StringIO
import csv

def main(argv):
  grep_term = ''
  seqfilename = ''
  try:
      opts, args = getopt.getopt(argv,"hg:s:",["seqfile=","grep="])
  except getopt.GetoptError:
    print 'Type RemoveSeqs.py -h for options'
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print 'RemoveSeqs.py -g <grep_term> -s <seqfile>'
       sys.exit()
    elif opt in ("-g", "--grep"):
       grep_term = arg
    elif opt in ("-s", "--seqfile"):
       seqfilename = arg
  out_seqs =[]
  for seq in  SeqIO.parse(seqfilename, "fasta"):
    if not grep_term in seq.id:
       out_seqs.append(seq)
  SeqIO.write(out_seqs, seqfilename, "fasta")
  
if __name__ == "__main__":
   main(sys.argv[1:])


