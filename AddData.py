#!/usr/local/bin/python

"""
A number of functions to add data to Selaginella database from csv files or from sequence
files"""


import sys, getopt, csv
import MySQLdb as mdb
from Bio import SeqIO
from os import path


def main(argv):
  usage = 'AddData -f <function> -i <inputfile>
  infilename = ''
  function = ''
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
        function? = arg

  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
    cur = con.cursor()
    if function == 'TAIR':
      add_TAIR(cur, infilename)
    if function == 'sequences':
      add_seqs(cur, infilename)
      
def add_TAIR(cur, infilename):
  infile = open(infilename, 'r')
  for line in infile.readlines():
    line = line.strip()
    fields = line.split('\t')
    locus = fields[0]
    if locus == 'Locus Identifier': #header row
      continue
    subfields = fields[2].split(';')
    description = subfields[0]
    try:
      symbol = fields[4]
    except IndexError:
      symbol = ''
    print "INSERT INTO tair (locus, description, symbol) VALUES ('%s', '%s', '%s')" % (locus, description, symbol)
    try:
      cur.execute("INSERT INTO tair (locus, description, symbol) VALUES (%s, %s, %s)", (locus, description, symbol))
    except mdb.IntegrityError, e:
      continue

def add_seqs(cur, infilename):
  species = path.split(infilename)[-1]
  species = species[:species.find('_')]
  print species
  for seq_record in SeqIO.parse(infilename, "fasta"):
    seq = seq_record.seq
    #id = seq_record.id
    id = "_".join(seq_record.id.split("|")[1:3])
    try:
      print 'INSERT INTO Sequences(seqid, sequence, species)  VALUES(%s, %s, %s)' % (id,seq,species,)
      cur.execute('INSERT INTO Sequences(seqid, sequence, species)  VALUES(%s, %s, %s)' , (id,seq,species,))
    except mdb.IntegrityError, e:
      warnings.warn("%s" % e)

if __name__ == "__main__":
   main(sys.argv[1:])


