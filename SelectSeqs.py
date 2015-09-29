#!/usr/local/bin/python

"""Select specified sequences from a larger file. This is useful, for example, when
   extracting a subclade from a larger phylogeny for more detailed analyses using a
   command like: nw_clade TREEFILE SEQ1 SEQ2 | nw_labels -I -"""
import sys, getopt, warnings

from Bio import SeqIO

def main(argv):
  usage = """SelectSeqs.py -i INFILE [-o OUTFILE] [-s SEQLIST]
             
             If OUTFILE not specified, sequences are written to STDOUT
             If SEQLIST not specified, sequences are read from SDTIN
          """
  infilename = ''
  outfilename = sys.stdout
  seqlistname = ''
  try:
      opts, args = getopt.getopt(argv,"hi:o:s:",["infile=", "outfile=", "seqlist="])
  except getopt.GetoptError:
    sys.exit(usage)
  for opt, arg in opts:
    if opt == "-h":
       sys.exit(usage)
    elif opt in ("-i", "--infile"):
       infilename = arg
    elif opt in ("-o", "--outfile"):
       outfilename = arg
    elif opt in ("-s", "--seqlist"):
       seqlistname = arg
  if not infilename:
     sys.exit(usage)
  
  if seqlistname:
     seqlist_fh = open(seqlistname, 'r')
  elif not sys.stdin.isatty():  #I added this so the script does not hang when seq list not provided at all
     seqlist_fh = sys.stdin
  else:
     sys.exit(usage)

  #make a dictionary of sequences to retain from seqlist
  seqlist = {}
  for line in seqlist_fh.readlines():
    seqlist[line.strip()] = 1
  
  #read through sequence file and add matching sequences to dictionary
  for seq in SeqIO.parse(infilename, "fasta"):
      if seq.id in seqlist:
         seqlist[seq.id] = seq
  
  #delete any dictionary entries that are not in the sequence file       
  for name, seq in seqlist.items():
     if seq == 1:
        warnings.warn("Sequence %s not in %s" % (name, infilename))
        del seqlist[name]
  
  #write selected sequences to outfile (or STDOUT)  
  SeqIO.write(seqlist.values(), outfilename, "fasta")
  
if __name__ == "__main__":
   main(sys.argv[1:])


