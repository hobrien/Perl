#!/usr/local/bin/python

import sys, getopt, csv
import MySQLdb as mdb
from os import path


def main(argv):
  infilename = ''
  function = ''
  usage = "GetData.py -f <function> -i <filename>"
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
      get_TAIR(cur, infilename)
      
def get_TAIR:
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
  

if __name__ == "__main__":
   main(sys.argv[1:])


