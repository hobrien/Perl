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
      fields = line.split(' ')
      cluster = fields[0]
      if cluster == 'Locus Identifier': #header row
        continue
      cluster = cluster.replace('cluster_', '')
      #print "SELECT seqid FROM sequences WHERE species = 'ATH' AND clusternum = %s" % (cluster)
      cur.execute("SELECT sequences.clusternum, sequences.seqid, tair.description, tair.symbol FROM sequences, tair WHERE sequences.seqid = tair.locus AND species = 'ATH' AND clusternum = %s", (cluster))
      for row in cur.fetchall():
        print '\t'.join(map(str, row))
      #cur.execute("INSERT INTO tair (locus, description, symbol) VALUES (%s, %s, %s)", (locus, description, symbol))
  

if __name__ == "__main__":
   main(sys.argv[1:])


