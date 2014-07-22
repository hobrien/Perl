#!/usr/local/bin/python

import sys, getopt, csv
import MySQLdb as mdb
from os import path
from Heathpy import flatten_GTF


def main(argv):
  infilename = ''
  outfilename = ''
  function = ''
  database = 'SelaginellaGenomics'
  usage = "GetData.py -f <function> -i <infilename> -s <species> -c <cluster> -d <database> -o <outfile>"
  try:
     opts, args = getopt.getopt(argv,"hf:i:s:c:d:o:",["function", "ifile=", "species=", "cluster=", "db=", "ofile="])
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
     elif opt in ("-s", "--species"):
        species = arg
     elif opt in ("-c", "--cluster"):
        cluster = arg
     elif opt in ("-d", "--db"):
        database = arg
     elif opt in ("-o", "--ofile"):
        outfilename = arg
  con = mdb.connect('localhost', 'root', '', database);

  if outfilename:
    outfile = open(outfilename, 'w')
  else : 
    outfile = sys.stdout
    
  with con:
    cur = con.cursor()
    if function == 'TAIR':
      get_TAIR(cur, infilename)
    elif function == 'sequences':
      get_seqs(cur, infilename)
    elif function == 'GTF':
      get_gtf(cur, infilename)
    elif function == 'seq_clusters':
      seq_clusters(cur, cluster, outfile)
    elif function == 'clusters':
      get_clusters(cur, species)
    elif function == 'locus':
      get_locus(cur, cluster)    

def get_locus(cur, locus):
  cur.execute("SELECT seqID, sequence FROM Sequences WHERE locusID = %s", (locus))
  for (accessionID, sequence) in cur.fetchall():
    print ">" + accessionID + "\n" + sequence

def get_clusters(cur, species):
  cur.execute("SELECT DISTINCT orthoID FROM OrthoGroups, CodingSequences, Sequences WHERE Orthogroups.geneID = CodingSequences.geneID AND CodingSequences.seqID = Sequences.seqID AND Sequences.species = %s", (species))
  for cluster in cur.fetchall():
    print cluster[0]
  

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
                  
def seq_clusters(cur, cluster, outfile):
  cur.execute("SELECT geneID FROM OrthoGroups WHERE orthoID = %s AND (geneID LIKE 'KRAUS%%' OR geneID LIKE 'MOEL%%' OR geneID LIKE 'UNC%%' OR geneID LIKE 'WILD%%')", cluster)
  if len(cur.fetchall()) > 0:
    cur.execute("SELECT DISTINCT CodingSequences.geneID, Sequences.sequence, CodingSequences.start, CodingSequences.end, CodingSequences.strand FROM Sequences, CodingSequences, OrthoGroups WHERE CodingSequences.geneID = OrthoGroups.geneID AND CodingSequences.seqID = Sequences.seqID AND OrthoGroups.orthoID = %s", (cluster))
    for (seqid, seq, start, end, strand) in cur.fetchall():
      seqid = seqid.replace("_", "|")
      outfile.write(">%s\n" % seqid)
      start = start - 1 #need to convert to python numbering
      if strand == '+':
        outfile.write(seq[start:end] + '\n')
      else:
        outfile.write(seq[start:end].reverse_complement() + '\n')

if __name__ == "__main__":
   main(sys.argv[1:])


