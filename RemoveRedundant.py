#!/usr/local/bin/python
"""

RemoveRedundant.py 14 Oct 2013

SYNOPSIS

RemoveRedundant.py -t <treefile> -s <seqfile> -n <species_name> -o <outfile> -l <logfile>

Options:
 -t treefile
 -n species name (KRAUS, UNC, WILD, MOEL)

DESCRIPTION
Takes a treefile and looks for monophyletic groups of sequences from a single species, 
then selects the sequence with the longest top blast hit and writes it to outfile

Consensus sequence is appended to the specified file and a summary of which sequences it
represents are appended to the specified logfile

NOTES

I've now got this set up to make a database for the blast hits so I can use this to easily
select the top hit and get the length.


AUTHOR

 Heath E. O'Brien (heath.obrien-at-gmail-dot-com)
"""

import sys
import getopt
import MySQLdb as mdb
import csv
from os import path
from ete2 import Tree  
from Bio import SeqIO
from Heathpy import parse_GTF
import glob

def main(argv):
  treefolder = ''
  gtf_filename = ''
  usage = 'RemoveRedundant.py -f <folder> -g <gtf_file>'
  try:
    opts, args = getopt.getopt(argv,"hf:g:",["folder=", "gtf="])
  except getopt.GetoptError:
    print usage
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print usage
      sys.exit()
    elif opt in ("-f", "--folder"):
      treefolder = arg
    elif opt in ("-g", "--gtf"):
      gtf_filename = arg

  #make dictionary to translate gene_ids back to contig names (this might be slow)
  #This will need to be a concatinated file if I'm going to loop over species
  if gtf_filename:
     gene_ids = get_id_dict(gtf_filename)

  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
    cur = con.cursor()
    for treefile in glob.glob(path.join(treefolder, '*')):
      if '.nwk' not in treefile:
        continue
      #read in tree and root by midpoint
      try:
        tree = Tree(treefile)
      except:
        continue
      root = tree.get_midpoint_outgroup()
      try:
        tree.set_outgroup(root)
      except:
        pass
      for species in ('KRAUS', 'MOEL', 'UNC', 'WILD'):
        remove_redundant(cur, tree, species, gene_ids)

def remove_redundant(cur, tree, species, gene_ids):
    species_leaves = []
    for leaf in tree:
        if species in leaf.name:
          species_leaves.append(leaf)
    
    for x in range(len(species_leaves)):
      for y in range(x+1, len(species_leaves)):
        if tree.get_distance(species_leaves[x], species_leaves[y], topology_only=True) <=2:
          seq1 = gene_ids[species_leaves[x].name]
          seq2 = gene_ids[species_leaves[y].name]
          if get_length(cur, seq1) >= get_length(cur, seq2):
            print "replaceing %s with %s" % (seq2, seq1)
            update_db(cur, seq2, seq1)
          else:
            print "replaceing %s with %s" % (seq1, seq2)
            update_db(cur, seq1, seq2)

"""This will find all sequences represented by a redundant sequence and update to the new rep seq, 
as well as updating the rep seq for the now-redundant seq"""
def update_db(cur, seqid, repseq):
  print seqid, repseq
  cur.execute("Select seqid from Sequences WHERE repseq = %s", (seqid))
  redundant_seqs =  cur.fetchall()
  for redundant_seq in redundant_seqs:
    print redundant_seq
    cur.execute("UPDATE Sequences SET repseq = %s WHERE seqid = %s", (repseq, redundant_seq[0]))
  cur.execute("UPDATE Sequences SET repseq = %s WHERE seqid = %s", (repseq, seqid))
  
def get_length(c, id):
  try:
    c.execute('SELECT qstart, qend FROM Blast_Selmo_all WHERE qseqid=%s', (id + '_0001',))
    return reduce(lambda x, y: y-x+1,c.fetchone())
  except TypeError:
    print "'SELECT qstart, qend FROM Blast_Selmo_all WHERE qseqid=%s': No Result" % id
    return 0

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
   
def get_id_dict(filename):
  seqids = {}
  for line in open(filename, 'r').readlines():
    feature = parse_GTF(line.split("\t"))
    if 'gene_id' in feature.keys() and feature['gene_id'] != 'NULL':
        seqids[feature['gene_id']] = feature['seqid']
  return seqids

def write_seqs(dict, file, format):
  fh = open(file, 'w')
  for seq in dict.values():
    fh.write(seq.format(format))
  fh.close()
  
if __name__ == "__main__":
   main(sys.argv[1:])
   #import timeit
   #print(timeit.timeit("main(sys.argv[1:])", setup="from __main__ import main", number=1))
