#!/usr/local/bin/python

import csv, sys, subprocess, argparse, os
from Bio import SeqIO

def parse_tblastn(args):
  assembly = args.seqfilename
  species_name = args.subject_name
  results = top_query(args)
  for result in results.values():
    subject = result[0]['sseqid']
    strand = result[0]['strand']
    header = ">" + species_name +"_" + result[0]['qseqid']
    seq = []
    for coord in combine_hits(result):
      seq.append(str(get_seq(assembly, args.index_filename, subject, coord[0], coord[1], strand)))
    print header
    print 'X'.join(seq)

def get_seq(seqfile, index, seqname, start, end, strand):
  sequence_db = SeqIO.index_db(index, seqfile, 'fasta')
  seq = sequence_db[seqname][start-1:end].seq
  if strand < 0:
    seq = seq.reverse_complement()
  seq = seq.translate()
  return seq  

def parse_blastp(argv):
  blastfilename = args.blastfilename   #List of top blast hits for each SUBJECT sequence
  assembly = args.seqfilename
  species_name = fullblastfilename.split('/')[-1].split('_')[0]
  evalue_cutoff = args.evalue
  with open(blastfilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader: 
      #Get blast stats and covert to numeric formats
      query, qlen, subject, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore = convert_blast_stats(row)

      if evalue > evalue_cutoff:
        continue
      header = ">Hygag_" + subject
      proc = subprocess.Popen(["GetSeq.pl %s %s" % (assembly, subject)], stdout=subprocess.PIPE, shell=True)
      (res, err) = proc.communicate()
      seq = ''.join(res.split('\n')[1:-1])
      print header
      print seq

def top_query(args):
  """returns a dictionary with subject sequences as keys and lists of all blast hsps for 
  the top-scoring query sequence as values
  -This means that, for example, if there are 2 related query sequences and homologous
  sequences on two contigs, the result will be (contig1:[hsp1, hsp2, hsp3], contig2:[hsps1, hsp2, hsp3])
  -Each value contains a list of all hsps from the top scoring query sequence for that contg
  -Note that this will exclude cases where multiple homologs are found on the same contig
  """  
  
  top_hits = {}
  blastfilename = args.blastfilename
  assembly = args.seqfilename
  species_name = args.subject_name
  evalue_cutoff = args.evalue
  with open(blastfilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader: 
      #Get blast stats and covert to numeric formats
      blast_result = parse_blast_stats(args.column_names, row)

      if blast_result['evalue'] > evalue_cutoff:
        continue
      if blast_result['sseqid'] not in top_hits:
        top_hits[blast_result['sseqid']] = [blast_result]
      elif blast_result['qseqid'] == top_hits[blast_result['sseqid']][0]['qseqid']:
        top_hits[blast_result['sseqid']].append(blast_result)
      elif blast_result['bitscore'] > top_hits[blast_result['sseqid']][0]['bitscore']:
        top_hits[blast_result['sseqid']] = [blast_result]
      else:
        pass 
  return top_hits

def parse_blast_stats(column_names, row):
  """This will convert stats from blast hits to the correct numeric format"""
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
  
def combine_hits(results):
  """This will return a sorted list of subject coordinate pairs with parts that overlap the
  query removed"""
  max_intron_length = 1000
  qstarts = []
  qends = []
  sstarts = []
  sends = []
  sframes = []
  coords = []
  strand = 1
  for result in results:
    if sstarts:
      if 'strand' * result['strand'] < 0: # checks if the strand is the same (neg * pos would be < 0)
        continue
    else:
      strand = result['strand']
    qstarts.append(result['qstart'])
    qends.append(result['qend'])
    sstarts.append(result['sstart'])
    sends.append(result['send'])
  #determine the correct sort order for hits
  sorted = qstarts[:] 
  sorted.sort()
  previous_qend = 0 #needed to figure out how much overlap to trim from subject sequences
  for qstart in sorted:
    index = qstarts.index(qstart) #this is the array postion for the current hit
    sstart = sstarts[index]
    send = sends[index]
    overlap = 0
    if qstarts[index] <= previous_qend:
      overlap =  (previous_qend - qstarts[index] + 1) * 3
    previous_qend = qends[index]
    if strand == 1:
      coords.append([sstart + overlap, send])
    else:
      coords.append([send, sstart - overlap])
  return coords
        
if __name__ == "__main__":
#=================== BEGIN ARGUMENT PARSING  =======================
  parser = argparse.ArgumentParser(description="Extract subject sequences from top scoring blast hit")
  parser.add_argument('blastfilename', help='name of blast outfile')
  parser.add_argument('seqfilename', help='Name of sequence file (fasta format)')
  parser.add_argument('--index_filename', '-i', dest='index_filename', default='',
                   help='Name of sequence index file (SQLLite3 database created by biopython)')
  parser.add_argument('--subject_name', '-n', dest='subject_name', default='',
                   help='Name of the subject organism')
  parser.add_argument('--evalue', '-e', dest='evalue', default=10, type=float,
                   help='maximum Evalue cutoff')
  parser.add_argument('--program', '-p', dest='program', default='', type=str,
                   help='blast program (tblastn or blastp)')  
  parser.add_argument('--outfmt', '-f', dest='column_names', default='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore', type=str,
                   help='blast output fields (default: standard -m 8 output)')  
  parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')
  args = parser.parse_args()
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
  if not args.program:
    if 'tblastn' in args.blastfilename:
      args.program = 'tblastn'
    elif 'tblastn' in args.blastfilename:
      args.program = 'blastp'
    else:
      print "can't determine blast algorithm from file name %s. Please specify with the --program option" % args.blastfilename
      print "type %s --help for more information" % os.path.basename(sys.argv[0])
      sys.exit()
  if not args.subject_name:
      args.subject_name = os.path.basename(args.seqfilename).split('.')[0]
  if not args.index_filename:
      args.index_filename = '.'.join(args.seqfilename.split('.')[:-1] + ['inx'])
      
#=================== END ARGUMENT PARSING  =======================
  if args.program =='tblastn':
    parse_tblastn(args)
  elif args.program == 'blastp':
    parse_blastp(args)
  else:
    print "blast program not recognized. Please specify either '--program blastp' or '--program tblastn'"
    print "type %s --help for more information" % os.path.basename(sys.argv[0])
    sys.exit()
