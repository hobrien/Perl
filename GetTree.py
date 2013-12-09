#!/usr/local/bin/python

import sys, getopt, csv
import MySQLdb as mdb
from os import path, system
from Heathpy import flatten_GTF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess

def main(argv):
  dirname = ''
  cluster = ''
  usage = "GetTree.py -c <cluster> -d <dirname>"
  try:
     opts, args = getopt.getopt(argv,"hc:d:",["cluster", "dir="])
  except getopt.GetoptError:
     print usage
     sys.exit(2)
  for opt, arg in opts:
     if opt == '-h':
        print usage
        sys.exit()
     elif opt in ("-c", "--cluser"):
        clusternum = int(arg)
     elif opt in ("-f", "--file"):
        dirname = arg
  tempfile = path.join(dirname, 'temp1')
  seq_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.fa')
  aln_file = path.join(dirname, 'translatorx_res.nt_ali.fasta')
  nexus_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.nex')
  tree_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.nwk')
  pdf_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.pdf')
  con = mdb.connect('localhost', 'root', '', 'Selaginella')
  with con:
    cur = con.cursor()
    seqs = get_seqs(cur, clusternum)
    SeqIO.write(seqs, seq_file, "fasta")
    print "translatorx_vLocal.pl -i %s -p F" % seq_file
    system("translatorx_vLocal.pl -i %s -p F" % seq_file)
    print "ConvertAln.py -i %s -o %s -f nexus" % (aln_file, nexus_file)
    system("ConvertAln.py -i %s -o %s -f nexus" % (aln_file, nexus_file))
    proc = subprocess.Popen(["trimal -gappyout -in %s -out %s -colnumbering" % (aln_file, tempfile)], stdout=subprocess.PIPE, shell=True)
    (trimal_res, err) = proc.communicate()
    add_exset(nexus_file, trimal_res)
    system("rm %s" % path.join(dirname, 'translatorx_*'))
    system("ConvertAln.py -i %s -f phylip -x fasta" % tempfile)
    system("rm %s" % tempfile)
    tempfile += '.phy'
    system("phyml -i %s" % tempfile)
    system("rm %s" % tempfile)
    system("mv %s %s" % (tempfile + '_phyml_tree.txt', tree_file))
    system("rm %s" % tempfile + '_phyml_stats.txt')
    system("ColourTree.py -t %s -o %s" % (tree_file, pdf_file))

def get_seqs(cur, clusternum):
  seqs = []
  cur.execute("SELECT seqid, sequence FROM Sequences WHERE repseq = 'NR' AND clusternum = %s", clusternum)
  for (seqid, sequence) in cur.fetchall():
    seq_record = SeqRecord(Seq(sequence), id=seqid)
    if seqid[:4] in ('KRAU', 'MOEL', 'UNCc', 'WILD'):  #BLUELEAF sequence. Need to remove non-coding portions
      cur.execute("SELECT start, end, strand, gene_id FROM orfs WHERE seqid = %s", seqid)
      for (start, end, strand, gene_id) in cur.fetchall():
        strand = int(strand)
        print seqid
        try:
          if int(gene_id.split('_')[1][1:]) != clusternum:
            print int(gene_id.split('_')[1][1:])
            raise ValueError
        except ValueError:
          print "skipping %s" % gene_id
          continue
        seq_record = seq_record[start-1:end]
        if strand == 1:
          name = '_'.join(map(str, [seqid, start, end]))
        elif strand == 0:
          name = '_'.join(map(str, [seqid, end, start]))
          seq_record = rev_com(seq_record)
        else:
          sys.exit('strand %s not recognized' % strand)
        seq_record.id = name 
    seq_record.description = ''
    seqs.append(seq_record)
  return seqs

def add_exset(file, trimal_res):
  import re
  from itertools import groupby, count
  incl = trimal_res.split(',')
  incl = map(int, incl)
  aln = open(file, 'r')
  for line in aln.readlines():
    match = re.search( r"nchar= *(\d+);", line, re.I)
    if match:
      length = int(match.group(1))
      break
  print length
  aln.close()
  excl = []
  for x in range(length):
    if not x in incl:
      excl.append(x+1)

  aln = open(file, 'a')
  aln.write("\nBEGIN ASSUMPTIONS;\n\t EXSET * TRIMAL  =  ")
  aln.write(' '.join(as_range(g) for _, g in groupby(excl, key=lambda n, c=count(): n-next(c))))  
  aln.write(";\n\nEND;\n")
  aln.write("BEGIN CODONS;\n\tCODONPOSSET * UNTITLED  =  1: 1 - %s\\3, 2: 2 - %s\\3, 3: 3 - %s\\3;\n\nEND;\n" % (length, length, length))
  aln.close
   
def as_range(iterable):
  l = list(iterable)
  if len(l) > 1:
    return '{0}-{1}'.format(l[0], l[-1])
  else:
    return '{0}'.format(l[0])

  

if __name__ == "__main__":
   main(sys.argv[1:])


