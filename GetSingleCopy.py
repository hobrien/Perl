#!/usr/local/bin/python

import sys
import MySQLdb as mdb

def main():
  con = mdb.connect('localhost', 'root', '', 'test');
  with con:
    cur = con.cursor()
    for clusternum in range(1,10000):
      for species in ('KRAUS', 'MOEL', 'UNC', 'WILD'):
        cur.execute("SELECT Sequences.seqid FROM Sequences, ClusterNum WHERE Sequences.seqid = ClusterNum.seqid AND Sequences.repseq = 'NR' AND Sequences.species = %s AND ClusterNum.cluster = %s", (species, clusternum))
        if len(cur.fetchall()) != 1:
          break
      else:
        print clusternum
  
if __name__ == "__main__":
   main()
