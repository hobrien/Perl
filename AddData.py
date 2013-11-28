#!/usr/local/bin/python

import sys, getopt, csv
import MySQLdb as mdb
from os import path


def main(argv):
  infilename = ''
  try:
     opts, args = getopt.getopt(argv,"hi:",["ifile="])
  except getopt.GetoptError:
     print 'AddOrthologs.py -i <inputfile>'
     sys.exit(2)
  for opt, arg in opts:
     if opt == '-h':
        print 'AddOrthologs.py -i <inputfile>'
        sys.exit()
     elif opt in ("-i", "--ifile"):
        infilename = arg

  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
    cur = con.cursor()
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


