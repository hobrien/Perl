#!/opt/local/bin/python

"""This is being rewritten to make a single exon feature for the length of the contig with a
single cds feature the length of the top blast hit. Chimeric transcripts are going to be
a big problem for this approach, but I can check for them independently (by comparing cds
and exon length) and break them up if necessary before doing the expression analysis
"""

import sys  
import getopt 
import csv
from Heathpy import flatten_GTF

def main(argv):
  blastfilename = ''
  gtf_filename = ''
  try:
    opts, args = getopt.getopt(argv,"hb:g:",["blast=","gtf="])
  except getopt.GetoptError:
    print 'Type Blast2GTF.py -h for options'
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print 'Blast2GTF.py -b <blastfile> -g <GTF_file>'
       sys.exit()
    elif opt in ("-g", "--gtf"):
       gtf_filename = arg
    elif opt in ("-b", "--blast"):
       blastfilename = arg

  previous_query = 'init'
  exon = {
  'source': 'Trinity',
  'feature': 'exon',
  'frame': '.',
  'score': '.'
  }
  cds = {
  'source': 'Ensembl',
  'feature': 'CDS',
  'frame': 1
  }
  with open(gtf_filename, 'wb') as outfile:
    gtf_writer = csv.writer(outfile, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
    with open(blastfilename, 'rU') as infile:
      reader=csv.reader(infile,delimiter='\t')
      for row in reader:
        qseqid, qlen, sacc, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore = row
        if pident < 70: #skip low-similarity hits
          continue

        #skip secondary hits from same query
        if qseqid == previous_query: #skip secondary hits
          continue 

        previous_query = qseqid
    
        #Get info for exon feature
        exon['seqname'] = qseqid
        exon['start'] = 1
        exon['end'] = qlen
        exon['gene_id'] =  sacc #This would be a very useful place to store cluster info
        exon['transcript_id'] =  exon['gene_id'] + '.1'
        if sstart > send:
          exon['strand'] = '+'
        else:
          exon['strand'] = '-'
        gtf_writer.writerow(flatten_GTF(exon))
  
        #get info for cds feature
        cds['seqname'] = exon['seqname']
        cds['score'] =  bitscore
        cds['gene_id'] =  exon['gene_id']
        cds['transcript_id'] =  exon['transcript_id']
        cds['strand'] = exon['strand']
        if cds['strand'] == '+':
          cds['start'] = qstart
          cds['end'] = qend
        else:
          cds['start'] = qend
          cds['end'] = qstart
        gtf_writer.writerow(flatten_GTF(cds))  

if __name__ == "__main__":
   main(sys.argv[1:])


