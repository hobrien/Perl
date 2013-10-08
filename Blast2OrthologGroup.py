#!/usr/local/bin/python
"""Use blast results to add sequences to ortholog groups
"""
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
    print 'Type Blast2OrthologGroups.py -h for options'
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print 'Blast2OrthologGroups.py -b <blastfile> -s <seqfile>'
       sys.exit()
    elif opt in ("-b", "--blastfile"):
       blastfilename = arg
    elif opt in ("-s", "--seqfile"):
       seqfilename = arg
  #Exclude sequences that cluster with chloroplast-encoded genes     
  chloroplast_clusters = ("Cluster_10175", "Cluster_10226", "Cluster_10232", "Cluster_10485", "Cluster_10984", "Cluster_11190", "Cluster_11826", "Cluster_11977", "Cluster_11980", "Cluster_12324", "Cluster_12326", "Cluster_12658", "Cluster_12697", "Cluster_1272", "Cluster_12984", "Cluster_13661", "Cluster_13666", "Cluster_13836", "Cluster_14312", "Cluster_14375", "Cluster_14565", "Cluster_1826", "Cluster_2867", "Cluster_29556", "Cluster_4023", "Cluster_4421", "Cluster_5258", "Cluster_5739", "Cluster_5972", "Cluster_615", "Cluster_7013", "Cluster_7860", "Cluster_8325", "Cluster_8384")     
  cluster_info = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "RefSeq", "SeqClusters.p")
  seq_groups = pickle.load( open( cluster_info, "rb" ) )
  indexfilename = ".".join(seqfilename.split(".")[:-1] + ["inx"])
  new_seqs = SeqIO.index_db(indexfilename, seqfilename, "fasta")
  used_seqs = {}
  with open(blastfilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader:
      qseqid, qlen, sacc, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore = row
      if sacc in seq_groups and not used_seqs.has_key(qseqid) and not seq_groups[sacc] in chloroplast_clusters:
        used_seqs[qseqid] = 1
        seq = new_seqs[qseqid]
        #cluster_filename = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "RefSeq", seq_groups[sacc] + ".fa")      
        cluster_filename = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "ContigClusters", seq_groups[sacc] + ".fa")     
        cluster_file = open(cluster_filename, "a")
        #print "saving %s to %s" % (qseqid, cluster_filename)
        id = seq.id
        if int(send) < int(sstart):
          seq = seq.reverse_complement()
          id = id + '_rc'
        seq.id, seq.description = id, id
        cluster_file.write(seq.format("fasta"))
        cluster_file.close()


if __name__ == "__main__":
   main(sys.argv[1:])


