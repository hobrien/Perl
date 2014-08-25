#!/usr/local/bin/python

import csv, sys, subprocess, argparse, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

"""
ReciprocalHit.py

25 Aug 2014

Heath O'Brien (heath.obrien-at-gmail-dot-com)

This produce a list of reciprocal best hits from two blast files (in tabular form)

This script requires biopython (see http://biopython.org/DIST/docs/install/Installation.html)

Type ReciprocalHit.py --help for instructions


"""

def reciprocal_hits(args):
  """returns a list of tuples containing reciprocal best hits"""
  reciprocal_hits = []
  results1 = top_hit(args.blastfilename1, args)
  results2 = top_hit(args.blastfilename2, args)
  for seqname in results1:
    seqmatch1 = results1[seqname]
    seqmatch2 = results2[seqmatch1]
    if seqmatch2 == seqname:
        reciprocal_hits.append((seqname, seqmatch1))
  return reciprocal_hits
      
def top_hit(blastfilename, args):
  """returns a dictionary with query sequences as keys and top hits as values
  """  
  
  top_hit_dict = {}
  evalue_cutoff = args.evalue
  with open(blastfilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader: 
      #Get blast stats and covert to numeric formats
      blast_result = parse_blast_stats(args.column_names, row)
      if blast_result['sseqid'] not in top_hit_dict:
        if blast_result['evalue'] <= evalue_cutoff:
          top_hit_dict[blast_result['sseqid']] = blast_result['qseqid']
  return top_hit_dict

def parse_blast_stats(column_names, row):
  """This will convert stats from blast hits to the correct numeric format and determine strand
  start and end coordinates are not reversed for negative strand results"""
  result = {}
  for field in column_names.split():
    try:
      if field in ('qlen', 'slen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qframe', 'sstart', 'send', 'sframe', 'sstrand'):
        result[field] = int(row.pop(0))
      elif field in ('pident', 'evalue', 'bitscore'):
        result[field] = float(row.pop(0))
      else:  
        result[field] = row.pop(0)
    except IndexError:
      sys.exit("number of columns does not match specified file format1. Please recheck column headers")
  if len(row) > 0:
    sys.exit("number of columns does not match specified file format. Please recheck column headers")
  if result['sstart'] < result['send']:
    result['strand'] = 1
  else:
    result['strand'] = -1
  
  return result    
   
def parse_args(input):
#=================== BEGIN ARGUMENT PARSING  =======================
  parser = argparse.ArgumentParser(description="Extract subject sequences from top scoring blast hit")
  parser.add_argument('blastfilename1', help='name of first blast outfile')
  parser.add_argument('blastfilename2', help='Name of second blast outfile')
  parser.add_argument('--evalue', '-e', dest='evalue', default=10, type=float,
                   help='maximum Evalue cutoff')
  parser.add_argument('--outfmt', '-f', dest='column_names', 
                   default='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', type=str,
                   help='blast output fields (default: standard -m 8 output)')  
  parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')
  args = parser.parse_args(input)
  args.column_names = args.column_names.replace('6 ', '')  #This is a hack that allows me to use the full custom fields specification from the the blast -outfmt command as input for this. 
  #substitute other query / subject identifiers for qseqid / sseqid if these are not present
  if not 'sseqid' in args.column_names:
    if 'sacc' in args.column_names:
       args.column_names = args.column_names.replace('sacc', 'sseqid')
    elif 'sgi' in args.column_names:
       args.column_names = args.column_names.replace('sgi', 'sseqid')
  if not 'qseqid' in args.column_names:
    if 'qacc' in args.column_names:
       args.column_names = args.column_names.replace('qacc', 'qseqid')
    elif 'qgi' in args.column_names:
       args.column_names = args.column_names.replace('qgi', 'qseqid')
  for field in ('sseqid', 'qseqid', 'bitscore', 'evalue', 'sstart', 'send', 'qstart', 'qend'):
    if field not in args.column_names:
      sys.exit("blast output must include %s" % field)
  args.blastfilename1 = os.path.abspath(args.blastfilename1)
  args.blastfilename2 = os.path.abspath(args.blastfilename2)
  
  return args
#=================== END ARGUMENT PARSING  =======================
 
        
if __name__ == "__main__":
  args = parse_args(sys.argv[1:])
  for reciprocal_hit in reciprocal_hits(args):
    print ' '.join(reciprocal_hit)
