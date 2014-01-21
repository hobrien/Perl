#!/usr/local/bin/python

import MySQLdb as mdb

def combine_clusters(cluster1, cluster2, db='Selaginella'):
  """replace the cluster number of cluster1 sequences with cluster2 (cluster1 sequences
     that are also in cluster2 are deleted"""
  con = mdb.connect('localhost', 'root', '', db)
  with con:
    cur = con.cursor()
    cur.execute("SELECT seqid FROM cluster_num WHERE cluster = %s", cluster1)
    for seq in cur.fetchall():
      cur.execute("SELECT seqid FROM cluster_num WHERE cluster = %s AND seqid = %s", (cluster2, seq[0]))
      if len(cur.fetchall()) > 0:
        #print "DELETE FROM cluster_num WHERE seqid = %s AND cluster = %s" % (seq[0], cluster1)
        cur.execute("DELETE FROM cluster_num WHERE seqid = %s AND cluster = %s", (seq[0], cluster1))
      else:
        cur.execute("UPDATE cluster_num SET cluster = %s WHERE seqid = %s", (cluster2, seq[0]))
