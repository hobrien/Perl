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
  logfilename = ''
  gtf_filename = ''
  nr_seq_filename = ''
  usage = 'RemoveRedundant.py -f <folder> -g <gtf_file> -l <logfile> -s <nr_seqs>'
  try:
    opts, args = getopt.getopt(argv,"hf:l:g:s:",["folder=","log=", "gtf=", "seqs="])
  except getopt.GetoptError:
    print usage
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print usage
      sys.exit()
    elif opt in ("-f", "--folder"):
      treefolder = arg
    elif opt in ("-l", "--log"):
      logfilename = arg
    elif opt in ("-g", "--gtf"):
      gtf_filename = arg
    elif opt in ("-s", "--seqs"):
      nr_seq_filename = arg

  #make dictionary to translate gene_ids back to contig names (this might be slow)
  #This will need to be a concatinated file if I'm going to loop over species
  if gtf_filename:
     gene_ids = get_id_dict(gtf_filename)

  #make dictionary of nr seqs (this is probably a temp
  if nr_seq_filename:
    nr_seqs = {}
    nr_seq_fh = open(nr_seq_filename, 'r')
    if gtf_filename:
      for seq in nr_seq_fh.readlines():
        nr_seqs[gene_ids[seq.rstrip()]] = 1
    else:
      for seq in nr_seq_fh.readlines():
        nr_seqs[seq] = 1
      
    nr_seq_fh.close()

  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
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

      #open corresponding fasta file and read sequences into dictionary
      seqfilename = treefile.replace("nwk", "fa")
      print seqfilename
      seqfilehandle = open(seqfilename, 'r')
      seq_dict = SeqIO.to_dict(SeqIO.parse(seqfilehandle, "fasta"))
      seqfilehandle.close()
      for key in seq_dict.keys():
        print key

      #cycle through monophyletic groups and select sequence with longest blast hit for each group
      for species_name in ('KRAUS', 'MOEL', 'UNC', 'WILD'):
        for node in t.get_monophyletic(values=[species_name], target_attr="species"):
          species_seqs = []
          for leaf in node:
            if 'kraussiana' not in leaf.name and 'willdenowii' not in leaf.name:  #exclude reference sequence from seq list
              species_seqs.append(leaf.name)
              seq_total += 1
          rep_seq = ''
          if nr_seq_filename:
            for seq in species_seqs:
              if seq in nr_seqs:
                rep_seq = seq
              else:
                del seq_dict[seq]
          else:
            for seq in species_seqs:
              if gtf_filename:
                try:
                  seq_name = gene_ids[seq]
                except KeyError:
                  sys.exit("%s in %s not found in GTF file" % (seq, treefile))
              else:
                seq_name = seq
              if rep_seq and get_length(cur, rep_seq) >= get_length(cur, seq_name): #saved rep seq longer than current seq
                del seq_dict[seq]
              else:                       #current seq longer (or no saved seq)
                if rep_seq:
                  del seq_dict[rep_seq]
                rep_seq = seq
          for seq in species_seqs:
            logfile.write("%s\n" % ", ".join([seq, rep_seq, cluster_num]))
      write_seqs(seq_dict, seqfilename, "fasta")
      print seq_total
    logfile.close()

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
        seqids[feature['seqid']] = feature['gene_id']
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
