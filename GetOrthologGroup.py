#!/opt/local/bin/python

import sys, getopt
import cPickle as pickle
from os import listdir, path, system
from Bio import SeqIO
from StringIO import StringIO
import csv

def main(argv):
  blastfilename = ''
  seqfilename = ''
  try:
      opts, args = getopt.getopt(argv,"hb:s:",["seqfile=","blastfile="])
  except getopt.GetoptError:
    print 'Type GetOrthologGroups.py -h for options'
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print 'GetOrthologGroups.py -b <blastfile> -s <seqfile>'
       sys.exit()
    elif opt in ("-b", "--blastfile"):
       blastfilename = arg
    elif opt in ("-s", "--seqfile"):
       seqfilename = arg
  species = path.split(seqfilename)[1].split("_")[0]
  print species
  cluster_info = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "RefSeq", "SeqClusters.p")
  seq_groups = pickle.load( open( cluster_info, "rb" ) )
  indexfilename = ".".join(seqfilename.split(".")[:-1] + ["inx"])
  new_seqs = SeqIO.index_db(indexfilename, seqfilename, "fasta")
  used_seqs = {}
  with open(blastfilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader:
      qseqid, qlen, sacc, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore = row
      if not used_seqs.has_key(qseqid):
        used_seqs[qseqid] = 1
        seq = new_seqs[qseqid]  
        #cluster_filename = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "RefSeq", seq_groups[sacc] + ".fa")      
        cluster_filename = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "RefSeq", "Homolog_groups", seq_groups[sacc] + ".fa")      
        cluster_file = open(cluster_filename, "a")
        #print "saving %s to %s" % (qseqid, cluster_filename)
        id = species + seq.id
        if int(send) < int(sstart):
          seq = seq.reverse_complement()
          id = id + '_rc'
        seq.id, seq.description = id, id
        cluster_file.write(seq.format("fasta"))
        cluster_file.close()


if __name__ == "__main__":
   main(sys.argv[1:])


