#!/usr/local/bin/python

import sys, subprocess
import MySQLdb as mdb


def main():
  con = mdb.connect('localhost', 'root', '', 'test');
  with con:
    cur = con.cursor()
    for clusternum in range(1,10000):
      sequences = []
      for species in ('KRAUS', 'MOEL', 'UNC', 'WILD'):
        cur.execute("SELECT Sequences.seqid, Sequences.sequence FROM Sequences, ClusterNum WHERE Sequences.seqid = ClusterNum.seqid AND Sequences.repseq = 'NR' AND Sequences.species = %s AND ClusterNum.cluster = %s", (species, clusternum))
        seqs = cur.fetchall()
        if len(seqs) != 1:
          break
        else:
          sequences.append(seqs)
      else:
        seqfile = open("seqfile.fa", 'w')
        for seq in sequences:
          seqfile.write(">%s\n%s\n" % seq[0])
        seqfile.close()
        aln_file = "Cluster_" + str(clusternum) + ".fa"
        subprocess.call(["mafft seqfile.fa >" + aln_file], shell=True)
    subprocess.call(["rm seqfile.fa"], shell=True)
if __name__ == "__main__":
   main()
