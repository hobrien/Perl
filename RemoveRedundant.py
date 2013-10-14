#!/usr/local/bin/python
"""

RemoveRedundant.py 14 Oct 2013

SYNOPSIS

RemoveRedundant.py -t <treefile> -s <seqfile> -n <species_name> -o <outfile> -l <logfile>

Options:
 -t treefile
 -s sequence file (fasta format)
 -n species name (KRAUS, UNC, WILD, MOEL)
 -o outfile for representative sequences (appended)
 -l logfile to record summary of which sequences are represented by sequences (appended)

DESCRIPTION
Takes a treefile and looks for monophyletic groups of sequences from a single species, 
then selects the sequence with the longest top blast hit and writes it to outfile

Consensus sequence is appended to the specified file and a summary of which sequences it
represents are appended to the specified logfile

NOTES

I've now got this set up to make a database for the blast hits so I can use this to easily
select the top hit and get the length.

I'm also going to index the seq file so I don't have to store that in memory either

AUTHOR

 Heath E. O'Brien (heath.obrien-at-gmail-dot-com)
"""

import sys
import getopt
import sqlite3
import csv
from os import path
from ete2 import Tree  
from Bio import SeqIO


def main(argv):
  treefile = ''
  seqfilename = ''
  outfilename = ''
  logfilename = ''
  species_name = ''
  blastfile = ''
  try:
    opts, args = getopt.getopt(argv,"ht:n:o:l:s:b:",["tree=","infile=","species=", "outfile=","log=","blast="])
  except getopt.GetoptError:
    print 'RemoveRedundant.py -t <treefile> -s <sequence_file> -n <species_name> -o <outfile> -l <logfile>'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print 'RemoveRedundant.py -t <treefile> -s <sequence_file> -n <species_name> -o <outfile> -l <logfile>'
      sys.exit()
    elif opt in ("-t", "--tree"):
      treefile = arg
    elif opt in ("-n", "--name"):
      species_name = arg
    elif opt in ("-o", "--outfile"):
      outfilename = arg
    elif opt in ("-l", "--log"):
      logfilename = arg
    elif opt in ("-s", "--sequences"):
      seqfilename = arg
    elif opt in ("-b", "--blast"):
      blastfile = arg

  #Build database of blast results (if one doesn't already exist)
  db_name = '.'.join(blastfile.split('.')[:-1] + ['db'])
  if not path.exists(db_name):
    make_db(blastfile, db_name)

  #Build database of sequences (if one doesn't already exist)   
  seq_db_name = '.'.join(seqfilename.split('.')[:-1] + ['inx'])
  seq_db = SeqIO.index_db(seq_db_name, seqfilename, "fasta")

  #read in tree and root by midpoint
  t = Tree(treefile)
  R = t.get_midpoint_outgroup()
  t.set_outgroup(R)

  #Add "species" feature to leaves on tree for easy lookup
  t = add_species(t)

  #parse tree file to get cluster number
  cluster_num = path.split(treefile)[1].split(".")[0] #get cluster number

  #cycle through monophyletic groups and select sequence with longest blast hit for each group
  conn = sqlite3.connect(db_name)
  c = conn.cursor()
  outfile = open(outfilename, "a")
  logfile = open(logfilename, "a")
  for node in t.get_monophyletic(values=[species_name], target_attr="species"):
    species_seqs = []
    for leaf in node:
      if 'kraussiana' not in leaf.name and 'willdenowii' not in leaf.name:  #exclude reference sequence from seq list
        species_seqs.append(leaf.name)
    rep_seq = ''
    for seq in species_seqs:
      if rep_seq and get_length(c, rep_seq.replace('_rc', '')) >= get_length(c, seq.replace('_rc', '')): #saved rep seq longer than current seq
        pass
      else:                       #current seq longer (or no saved seq)
        rep_seq = seq
    if rep_seq:
      rep_seq = 'KRAUScomp118232_c0_seq1_rc'
      if rep_seq.find('_rc') == -1:
        outfile.write(seq_db.get_raw(rep_seq))
      else:
        outseq = seq_db[rep_seq.replace('_rc','')].reverse_complement()
        outseq.id = rep_seq
        outseq.description = rep_seq
        outfile.write(outseq.format("fasta"))
      for seq in species_seqs:
        logfile.write("%s\n" % ", ".join([seq, rep_seq, cluster_num]))

  logfile.close()
  outfile.close()

 
def make_db(blastfile, db_name):
  conn = sqlite3.connect(db_name)
  c = conn.cursor()
  #strand: 1 = pos, 0 = neg
  c.execute('''CREATE TABLE hits
                (qseqid text, qlen int, sacc text, slen int, pident read, length int, mismatch int, gapopen int, qstart int, qend int, qframe int, sstart int, send int, sframe int, evalue float, bitscore float, strand bit, hitnum int)''')
  with open(blastfile, 'rU') as infile:
    reader=csv.reader(infile,delimiter='\t')
    for row in reader:
      row.append(1)  #add strand info to row
      if int(row[8]) > int(row[9]):               #qstart > qend
        (row[8], row[9]) = (row[9], row[8])       #reverse qstart and qend
        row[-1] = 0                               #change strand to neg
      if int(row[11]) > int(row[12]):             #sstart > send
        (row[11], row[12]) = (row[12], row[11])   #reverse sstart and send
        row[-1] = 0                               #change strand to neg
      c.execute('SELECT COUNT(*) FROM hits WHERE qseqid=?', (row[0],)) #count number of entries for query in db
      row.append(c.fetchone()[0])                                      #add count num to row
      c.execute("INSERT INTO hits VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,?)", row)
    conn.commit()
  conn.close()

def get_length(c, id):
  c.execute('SELECT qstart, qend FROM hits WHERE qseqid=? and hitnum=0', (id,))
  return reduce(lambda x, y: y-x+1,c.fetchone())

def add_species(tree): 
   species_list = ['UNC', 'MOEL', 'WILD', 'KRAUS']
   for leaf in tree:
     to_add = "other"
     for species in species_list:
       if species in leaf.name:
         to_add = species
     if "kraussiana" in leaf.name:       #Include reference sequences from 1kp project
       to_add = "KRAUS"
     elif "willdenowii" in leaf.name:
       to_add = "WILD"
     leaf.add_features(species=to_add)
   return tree
   

if __name__ == "__main__":
   main(sys.argv[1:])


