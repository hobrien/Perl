#!/usr/local/bin/python

import sys, getopt, csv
import MySQLdb as mdb
from os import path
from Heathpy import flatten_GTF
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


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
    
  with con:
    cur = con.cursor()
    if function == 'TAIR':
      get_TAIR(cur, infilename)
    elif function == 'sequences':
      get_seqs(cur, infilename)
    elif function == 'GTF':
      get_gtf(cur, species)
    elif function == 'seq_clusters':
      seq_clusters(cur, cluster, outfilename)
    elif function == 'nr_clusters':
      nr_clusters(cur, cluster, outfilename)
    elif function == 'clusters':
      get_clusters(cur, species)
    elif function == 'locus':
      get_locus(cur, cluster)    
    elif function == 'nr_seqs':
      get_nr_seqs(cur, species)    
    elif function == 'nr':
      get_nr(cur, species)    
    elif function == 'counts':
      get_counts(cur, species)    

def get_counts(cur, species):
  if species == 'all':
    print '\t'.join(('cluster', 'KRAUS1', 'KRAUS2', 'KRAUS3', 'KRAUS4', 'MOEL1', 'MOEL2', 'MOEL3', 'UNC1', 'UNC2', 'UNC3', 'UNC4', 'WILD1', 'WILD2', 'WILD3', 'WILD4'))
  elif species == 'MOEL':
    print '\t'.join(('cluster', species + '1',  species + '2',  species + '3'))
  else:
    print '\t'.join(('cluster', species + '1',  species + '2',  species + '3',  species + '4'))
  cur.execute("SELECT DISTINCT orthoID FROM OrthoGroups WHERE non_redundant = 1")
  for row in cur.fetchall():
    counts = ['cluster_' + str(row[0])]
    if species == 'all':
      cur.execute("SELECT Counts.* FROM Counts, OrthoGroups WHERE Counts.geneID = OrthoGroups.geneID AND OrthoGroups.orthoID = %s AND OrthoGroups.non_redundant = 1", (row[0]))
    else:
      cur.execute("SELECT Counts.* FROM Counts, OrthoGroups WHERE Counts.geneID = OrthoGroups.geneID AND OrthoGroups.orthoID = %s AND OrthoGroups.non_redundant = 1 AND Counts.geneID LIKE %s", (row[0], species + '%')) 
    for results in cur.fetchall():
      counts.append(int(results[1])+int(results[5]))
      counts.append(int(results[2])+int(results[6]))
      counts.append(int(results[3])+int(results[7]))
      if not 'MOEL' in results[0]:
        counts.append(int(results[4])+int(results[8]))
    if len(counts) == 16 or (species != 'all' and len(counts) > 1):
      print '\t'.join(map(str, counts))
    
def get_nr_seqs(cur, species):
  cur.execute("SELECT Orthogroups.orthoID, Sequences.seqID, Sequences.sequence FROM Sequences, CodingSequences, OrthoGroups WHERE Sequences.seqID = CodingSequences.seqID AND CodingSequences.geneID = OrthoGroups.geneID AND OrthoGroups.non_redundant=1 AND Sequences.species = %s", (species))
  for (cluster, seqID, sequence) in cur.fetchall():
    print ">" + seqID + "\n" + sequence

def get_nr(cur, species):
  cur.execute("SELECT Orthogroups.orthoID, Sequences.seqID, Sequences.sequence FROM Sequences, CodingSequences, OrthoGroups WHERE Sequences.seqID = CodingSequences.seqID AND CodingSequences.geneID = OrthoGroups.geneID AND OrthoGroups.non_redundant=1 AND Sequences.species = %s", (species))
  for (cluster, seqID, sequence) in cur.fetchall():
    print "cluster_" + str(cluster)

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
    
def get_gtf(cur, species):
  #print "SELECT seqID, start_pos, end_pos, strand, geneID FROM CodingSequences WHERE species = %s" % (species)
  cur.execute("SELECT seqID, start_pos, end_pos, strand, geneID FROM CodingSequences WHERE species = %s", (species))
  for row in cur.fetchall():
      gtf = { 'seqname': row[0],
              'source': 'transdecoder',
              'feature': 'CDS',
              'start': row[1],
              'end': row[2],
              'strand': row[3],
              'frame': '.',
              'gene_id': row[4],
              'transcript_id': row[4] + '.1'
            }
      print flatten_GTF(gtf)  
                  
def seq_clusters(cur, cluster, outfilename):
  cur.execute("SELECT geneID FROM OrthoGroups WHERE orthoID = %s AND (geneID LIKE 'KRAUS%%' OR geneID LIKE 'MOEL%%' OR geneID LIKE 'UNC%%' OR geneID LIKE 'WILD%%')", cluster)
  if len(cur.fetchall()) > 0:
    if outfilename:
      outfile = open(outfilename, 'w')
    else : 
      outfile = sys.stdout
    cur.execute("SELECT DISTINCT CodingSequences.geneID, Sequences.sequence, CodingSequences.start, CodingSequences.end, CodingSequences.strand FROM Sequences, CodingSequences, OrthoGroups WHERE CodingSequences.geneID = OrthoGroups.geneID AND CodingSequences.seqID = Sequences.seqID AND OrthoGroups.orthoID = %s", (cluster))
    for (seqid, seq, start, end, strand) in cur.fetchall():
      seqid = seqid.replace("_", "|")
      seq = Seq(seq, IUPAC.unambiguous_dna)
      outfile.write(">%s\n" % seqid)
      #print ">%s" % seqid
      start = start - 1 #need to convert to python numbering
      if strand == '+':
        outfile.write(str(seq[start:end]) + '\n')
        #print seq[start:end]
      else:
        outfile.write(str(seq[start:end].reverse_complement()) + '\n')
        #print seq[start:end].reverse_complement()

def nr_clusters(cur, cluster, outfilename):
  cur.execute("SELECT geneID FROM OrthoGroups WHERE orthoID = %s AND non_redundant = 1 AND (geneID LIKE 'KRAUS%%' OR geneID LIKE 'MOEL%%' OR geneID LIKE 'UNC%%' OR geneID LIKE 'WILD%%')", cluster)
  if len(cur.fetchall()) > 0:
    if outfilename:
      outfile = open(outfilename, 'w')
    else : 
      outfile = sys.stdout
    cur.execute("SELECT DISTINCT CodingSequences.geneID, Sequences.sequence, CodingSequences.start, CodingSequences.end, CodingSequences.strand FROM Sequences, CodingSequences, OrthoGroups WHERE CodingSequences.geneID = OrthoGroups.geneID AND CodingSequences.seqID = Sequences.seqID AND non_redundant = 1 AND  OrthoGroups.orthoID = %s", (cluster))
    for (seqid, seq, start, end, strand) in cur.fetchall():
      seqid = seqid.replace("_", "|")
      seq = Seq(seq, IUPAC.unambiguous_dna)
      outfile.write(">%s\n" % seqid)
      #print ">%s" % seqid
      start = start - 1 #need to convert to python numbering
      if strand == '+':
        outfile.write(str(seq[start:end]) + '\n')
        #print seq[start:end]
      else:
        outfile.write(str(seq[start:end].reverse_complement()) + '\n')
        #print seq[start:end].reverse_complement()

if __name__ == "__main__":
   main(sys.argv[1:])


