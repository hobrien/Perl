#!/usr/local/bin/python

import sys, getopt, csv
import MySQLdb as mdb
from os import path
from Heathpy import flatten_GTF


def main(argv):
  infilename = ''
  function = ''
  usage = "GetData.py -f <function> -i <filename>"
  try:
     opts, args = getopt.getopt(argv,"hf:i:",["function", "ifile="])
  except getopt.GetoptError:
     print usage
     sys.exit(2)
  for opt, arg in opts:
     if opt == '-h':
        print usage
        sys.exit()
     elif opt in ("-i", "--ifile"):
        infilename = arg
     elif opt in ("-f", "--function"):
        function = arg

  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
    cur = con.cursor()
    if function == 'TAIR':
      get_TAIR(cur, infilename)
    elif function == 'sequences':
      get_seqs(cur, infilename)
    elif function == 'GTF':
      get_gtf(cur, infilename)
    elif function == 'seq_clusters':
      seq_clusters(cur, infilename)
      
def get_TAIR(cur, infilename):
  infile = open(infilename, 'r')
  for line in infile.readlines():
    line = line.strip()
    fields = line.split(' ')
    cluster = fields[0]
    if cluster == 'Locus Identifier': #header row
      continue
    cluster = cluster.replace('cluster_', '')
    #print "SELECT seqid FROM sequences WHERE species = 'ATH' AND clusternum = %s" % (cluster)
    cur.execute("SELECT sequences.clusternum, sequences.seqid, tair.description, tair.symbol FROM sequences, tair WHERE sequences.seqid = tair.locus AND species = 'ATH' AND clusternum = %s", (cluster))
    for row in cur.fetchall():
      print '\t'.join(map(str, row))
  
def get_seqs(cur, name):
  cur.execute("SELECT seqid, sequence FROM Sequences WHERE repseq = 'NR' AND species = %s", (name))
  for (seqid, seq) in cur.fetchall():
    print ">%s" % seqid
    print seq
    
def get_gtf(cur, name):
  cur.execute("SELECT seqid, sequence FROM Sequences WHERE repseq = 'NR' AND species = %s", (species))
  for (seqid, seq) in cur.fetchall():
    cur.execute("SELECT * FROM orfs WHERE seqid = %s", (seqid,))
    for row in cur2.fetchall():
      gtf = { 'seqname': row[1],
              'source': row[2],
              'feature': row[3],
              'start': row[4],
              'end': row[5],
              'strand': row[6],
              'frame': row[7],
              'score': row[8],
              'gene_id': row[9],
              'transcript_id': row[10]
            }
      if row[11]:
        gtf['note'] = row[11]
      print flatten_GTF(gtf)  
                  
def seq_clusters(cur, name):
  if name == 'ATH':
    cur.execute("SELECT seqid, sequence, clusternum FROM Sequences WHERE species = 'ATH'")
  else:
    cur.execute("SELECT seqid, sequence, clusternum FROM Sequences WHERE repseq = 'NR' AND species = %s", (name))
  for (seqid, seq, clusternum) in cur.fetchall():
    print clusternum
    if clusternum == None:
      continue
    cluster_filename = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "NR_Clusters", "Cluster_" + str(clusternum) + ".fa")
    cluster_file = open(cluster_filename, "a")
    cluster_file.write(">%s\n%s\n" % (seqid, seq))
    cluster_file.close()
        

if __name__ == "__main__":
   main(sys.argv[1:])


