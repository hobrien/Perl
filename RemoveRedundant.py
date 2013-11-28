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

import sys, warnings
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


  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
    cur = con.cursor()
    for treefile in glob.glob(path.join(treefolder, '*')):
      if '.nwk' not in treefile:
        continue
      #read cluster number from filename
      cluster = int(treefile.split('_')[1].split('.')[0])
      print cluster
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
      tree = add_species(tree)
      for leaf in tree:
        if 'scaffold' in leaf.name:
          seq = "_".join(leaf.name.split("-")[1:3])
          try:
            cur.execute("SELECT clusternum FROM Sequences WHERE seqid = %s", ( seq ))
            saved_cluster = int(cur.fetchone()[0])
          except:
            saved_cluster = cluster
          if saved_cluster != cluster:
            warnings.warn("Sequence %s is in the cluster %s tree, but listed as cluster %s in the database" % (seq, cluster, saved_cluster))
          cur.execute("UPDATE Sequences SET clusternum = %s WHERE seqid = %s", ( cluster, seq ))
      for species in ('KRAUS', 'MOEL', 'UNC', 'WILD'):
        remove_redundant(cur, tree, species)

def remove_redundant(cur, tree, species):
    species_clades = []
    for node in tree.get_monophyletic(values=[species], target_attr="species"):
      species_clades.append(node)

    #Allow at most one clade of sequences from other species
    redundant = {}   
    for x in range(len(species_clades)):
      sequences = []
      for leaf in species_clades[x]:
        sequences.append(leaf.name)
      for y in range(x+1, len(species_clades)):
        if tree.get_distance(species_clades[x], species_clades[y], topology_only=True) <=2:
          for leaf in species_clades[y]:
            sequences.append(leaf.name)
      #print "Finding NR sequences for %s" % ', '.join(sequences)
      nr_seq = ''
      for gene_id in sequences:
        if 'scaffold' in gene_id:
          seq = "_".join(gene_id.split("-")[1:3])
        else:
          cur.execute("SELECT Sequences.seqid, Sequences.repseq FROM Sequences, orfs WHERE orfs.seqid=Sequences.seqid AND orfs.gene_id = %s", (gene_id))
          try:
            (seq, repseq) = cur.fetchone()
          except TypeError:
            sys.exit("no result for %s" % gene_id)
          if repseq != 'NR':
            print "Sequence %s repreesented by %s in the database but present in the tree" % (gene_id, repseq)
            pass
        if nr_seq:
          if get_length(cur, seq) > get_length(cur, nr_seq):
            update_db(cur, nr_seq, seq)
            nr_seq = seq
          else:
            update_db(cur, seq, nr_seq)
        else:
          nr_seq = seq

"""This will find all sequences represented by a redundant sequence and update to the new rep seq, 
as well as updating the rep seq for the now-redundant seq"""
def update_db(cur, seqid, repseq):
  #Get cluster number and species info from db because these columns are indexed for fast retrieval
  #This assumes that sequences in different clusters are not redundant to each other, which seems safe
  cur.execute("SELECT clusternum, species FROM Sequences WHERE seqid = %s", (seqid))
  (clusternum, species) = cur.fetchone()
  
  cur.execute("SELECT seqid FROM Sequences WHERE repseq = %s AND clusternum = %s AND species = %s", (seqid, clusternum, species))
  redundant_seqs =  cur.fetchall()
  for redundant_seq in redundant_seqs:
    cur.execute("UPDATE Sequences SET repseq = %s WHERE seqid = %s", (repseq, redundant_seq[0]))
  cur.execute("UPDATE Sequences SET repseq = %s WHERE seqid = %s", (repseq, seqid))
  cur.execute("UPDATE Sequences SET repseq = 'NR' WHERE seqid = %s", (repseq))

def get_length(c, id):
  try:
    c.execute('SELECT qstart, qend FROM Blast_Selmo_all WHERE qseqid=%s', (id + '_0001',))
    return reduce(lambda x, y: y-x+1,c.fetchone())
  except TypeError:
    #print "'SELECT qstart, qend FROM Blast_Selmo_all WHERE qseqid=%s': No Result" % id
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
     
if __name__ == "__main__":
   main(sys.argv[1:])
   #import timeit
   #print(timeit.timeit("main(sys.argv[1:])", setup="from __main__ import main", number=1))
