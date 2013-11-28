#!/usr/local/bin/python

"""
A number of functions to add data to Selaginella database from csv files or from sequence
files"""


import sys, getopt, csv
import MySQLdb as mdb
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
      AddTAIR(cur, infilename)
      
def addTAIR(cur, infilename):
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

if __name__ == "__main__":
   main(sys.argv[1:])


