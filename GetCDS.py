#!/usr/local/bin/python

import sys, getopt
from Bio import Entrez, SeqIO
Entrez.email = "heath.obrien@gmail.com"


def main(argv):
   inputfile = ''
#   outputfile = ''
   usage = 'GetCDS.py -i <inputfile>'
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print usage
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print usage
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
   for pep in SeqIO.parse(inputfile, "genbank"):
      cds = get_cds(pep)
      print cds.format("fasta")

def get_cds(pep):  
   strand = 1
   seq = ''
   coords = ''
   for feature in pep.features:
      if feature.type == "CDS":
         coords = feature.qualifiers['coded_by'][0]
         #tag = feature.qualifiers['locus_tag'][0]
   if coords.find('complement') == 0:
      coords = coords[11:-1]
      strand = 2
   if coords.find('join') == 0:
      coords = coords[5:-1]
   coords = coords.replace('>', '')
   coords = coords.replace('<', '')
   for interval in coords.split(', '):
      id = interval.split(':')[0]
      start = interval.split(':')[1].split('..')[0]
      end = interval.split(':')[1].split('..')[1]
      handle = Entrez.efetch(db="nucleotide", 
                       id=id, 
                       rettype="fasta", 
                       strand=strand, 
                       seq_start=start, 
                       seq_stop=end)
      record = SeqIO.read(handle, "fasta")
      handle.close()
      if strand == 1:
         seq = seq + record.seq
      else:
         seq = record.seq + seq
   pep.seq = seq
   #pep.id = tag
   return pep

if __name__ == "__main__":
   main(sys.argv[1:])


