#!/usr/local/bin/python

import sys
import MySQLdb as mdb
from string import maketrans
from os import path

def main():
  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
    used_names = {}
    cur = con.cursor()
    cur.execute("SELECT Sequences.seqid, Sequences.sequence, Sequences.clusternum, orfs.start, orfs.end, orfs.strand, orfs.gene_id FROM Sequences, orfs WHERE Sequences.seqid = orfs.seqid AND Sequences.repseq='NR'")
    for (seqid, seq, clusternum, start, end, strand, gene_id) in cur.fetchall():
#      clusternum = int(clusternum)
#      start = int(start)
#      end = int(end)
      strand = int(strand)
      print seqid
      try:
        if int(gene_id.split('_')[1][1:]) != clusternum:
          raise ValueError
      except ValueError:
        print "skipping %s" % gene_id
        continue
     
      seq = seq[start-1:end]
      if strand == 1:
        name = '_'.join(map(str, [seqid, start, end]))
      elif strand == 0:
        name = '_'.join(map(str, [seqid, end, start]))
        seq = rev_com(seq)
      else:
        sys.exit('strand %s not recognized' % strand)   
      cluster_filename = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "NR_Clusters", "Cluster_" + str(clusternum) + ".fa")
      cluster_file = open(cluster_filename, "a")
      cluster_file.write(">%s\n%s\n" % (name, seq))
      cluster_file.close()

def rev_com(seq):
  seq = seq.upper()
  seq = seq.translate(maketrans('ACGT', 'TGCA'))
  seq = seq[::-1]
  return seq
  
if __name__ == "__main__":
  main()