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
  seqfilename = ''
  outfilename = ''
  logfilename = ''
  species_name = ''
  gtf_filename = ''
  
  try:
    opts, args = getopt.getopt(argv,"ht:n:o:l:s:b:g:",["tree=","infile=","species=", "outfile=","log=", "gtf="])
  except getopt.GetoptError:
    print 'RemoveRedundant.py -t <treefile> -s <sequence_file> -n <species_name> -o <outfile> -l <logfile>'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print 'RemoveRedundant.py -t <treefile> -s <sequence_file> -n <species_name> -o <outfile> -l <logfile>'
      sys.exit()
    elif opt in ("-t", "--tree"):
      treefolder = arg
    elif opt in ("-n", "--name"):
      species_name = arg
    elif opt in ("-o", "--outfile"):
      outfilename = arg
    elif opt in ("-l", "--log"):
      logfilename = arg
    elif opt in ("-s", "--sequences"):
      seqfilename = arg
    elif opt in ("-g", "--gtf"):
      gtf_filename = arg

  #Build database of sequences (if one doesn't already exist)   
  indexfilename = seqfilename + ".inx"
  seq_db = SeqIO.index_db(indexfilename, seqfilename, "fasta")
  
  #make dictionary to translate gene_ids back to contig names (this might be slow)
  if gtf_filename:
     gene_ids = get_id_dict(gtf_filename)


  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
    outfile = open(outfilename, "w")
    logfile = open(logfilename, "w")
    cur = con.cursor()
    seq_total = 0
    for treefile in glob.glob(path.join(treefolder, '*')):
      if '.nwk' not in treefile:
        continue 
      #read in tree and root by midpoint
      try:
        t = Tree(treefile)
      except:
        continue
      R = t.get_midpoint_outgroup()
      try:
        t.set_outgroup(R)
      except:
        pass
      #Add "species" feature to leaves on tree for easy lookup
      t = add_species(t)

      #parse tree file to get cluster number
      cluster_num = path.split(treefile)[1].split(".")[0] #get cluster number

      #cycle through monophyletic groups and select sequence with longest blast hit for each group
      for node in t.get_monophyletic(values=[species_name], target_attr="species"):
        species_seqs = []
        for leaf in node:
          if 'kraussiana' not in leaf.name and 'willdenowii' not in leaf.name:  #exclude reference sequence from seq list
            if gtf_filename:
              try:
                species_seqs.append(gene_ids[leaf.name])
              except KeyError:
                print "%s in %s not found in GTF file" % (leaf.name, treefile)
            else:
              species_seqs.append(leaf.name)
            seq_total += 1
        rep_seq = ''
        for seq in species_seqs:
          if rep_seq and get_length(cur, rep_seq) >= get_length(cur, seq): #saved rep seq longer than current seq
            pass
          else:                       #current seq longer (or no saved seq)
            rep_seq = seq
        if rep_seq:
          try:
            outfile.write(seq_db.get_raw(rep_seq))
          except KeyError:
            sys.exit("No sequence %s in %s" % ( rep_seq, seqfilename))
          for seq in species_seqs:
            logfile.write("%s\n" % ", ".join([seq, rep_seq, cluster_num]))


    print seq_total
    logfile.close()
    outfile.close()

def get_length(c, id):
  try:
    c.execute('SELECT qstart, qend FROM Blast WHERE qseqid=%s and hitnum=1', (id,))
    return reduce(lambda x, y: y-x+1,c.fetchone())
  except TypeError:
    print "'SELECT qstart, qend FROM Blast WHERE qseqid=%s and hitnum=1': No Result" % id
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

if __name__ == "__main__":
   main(sys.argv[1:])
   #import timeit
   #print(timeit.timeit("main(sys.argv[1:])", setup="from __main__ import main", number=1))
