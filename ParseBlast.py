#!/usr/local/bin/python
"""This script"""

import csv, sys, subprocess, argparse, os

def parse_tblastn(argv):
  fullblastfilename = argv[0]   #List of top blast hits for each SUBJECT sequence
  uniqueblastfilename = argv[1]
  assembly = argv[2]
  species_name = fullblastfilename.split('/')[-1].split('_')[0]
  evalue_cutoff = float(argv[3])
  with open(uniqueblastfilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader: 
      #Get blast stats and covert to numeric formats
      query, qlen, subject, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore = convert_blast_stats(row)

      if evalue > evalue_cutoff:
        continue
      #Get a list of all hits for the same query and subject (this is a bit of a hack, but it works)
      proc = subprocess.Popen(["grep '%s' %s | grep %s" % (query, fullblastfilename, subject)], stdout=subprocess.PIPE, shell=True)
      (blast_res, err) = proc.communicate()
      header = ">" + species_name +"_" + query
      seq = []
      for coord in combine_hits(blast_res.split('\n')[:-1]):
        if sframe < 0:
          proc = subprocess.Popen(["GetSeq.pl %s %s | trunc.pl %s %s | rc.pl | translate.pl" % tuple([assembly, subject] + coord)], stdout=subprocess.PIPE, shell=True)
        else:
          proc = subprocess.Popen(["GetSeq.pl %s %s | trunc.pl %s %s | translate.pl" % tuple([assembly, subject] + coord)], stdout=subprocess.PIPE, shell=True)
        (res, err) = proc.communicate()
        seq.append(''.join(res.split('\n')[1:-1]))
      print header
      print 'X'.join(seq)

def parse_blastp(argv):
  uniqueblastfilename = argv[0]   #List of top blast hits for each SUBJECT sequence
  assembly = argv[1]
  
  evalue_cutoff = float(argv[2])
  with open(uniqueblastfilename, 'rU') as f:
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

def convert_blast_stats(result):
  """This will convert stats from blast hits to the correct numeric formate"""
  query, qlen, subject, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore = result
  qlen, slen, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe = map(lambda x: int(x), (qlen, slen, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe))
  pident, evalue, bitscore = map(lambda x: float(x), (pident, evalue, bitscore))
  return query, qlen, subject, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore
  
def combine_hits(hits):
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
  for line in hits:
    #Get blast stats and covert to numeric formats
    query, qlen, subject, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore = convert_blast_stats(line.split('\t'))
    #make lists of coordinates for query and subject sequences
    if sstarts:
      if strand * sframe < 0: #the second part of this checks if the strand is the same (neg * pos would be < 0)
        continue
    if sframe < 0: strand = -1
    qstarts.append(qstart)
    qends.append(qend)
    sstarts.append(sstart)
    sends.append(send)
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
  parser = argparse.ArgumentParser(description="Extract subject sequences from top scoring blast hit")
  parser.add_argument('filename', help='name of blast outfile')
  parser.add_argument('--evalue', '-e', dest='evalue', default=10, type=float,
                   help='maximum Evalue cutoff')
  parser.add_argument('--program', '-p', dest='program', default='', type=str,
                   help='blast program (tblastn or blastp)')  
  parser.add_argument('--out_format', '-f', dest='column_names', default='qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore', type=str,
                   help='blast output fields (default: standard -m 8 output)')  
  parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')
  args = parser.parse_args()
  if not args.program:
    if 'tblastn' in args.filename:
      args.program = 'tblastn'
    elif 'tblastn' in args.filename:
      args.program = 'blastp'
    else:
      print "can't determine blast algorithm from file name %s. Please specify with the --program option" % args.filename
      print "type %s --help for more information" % os.path.basename(sys.argv[0])
      sys.exit()
  if args.program =='tblastn':
    parse_tblastn(args)
  elif args.program == 'blastp':
    parse_blastp(args)
  else:
    print "blast program not recognized. Please specify either '--program blastp' or '--program tblastn'"
    print "type %s --help for more information" % os.path.basename(sys.argv[0])
    sys.exit()
