#!/opt/local/bin/python

"""This is being rewritten to make a single cds feature the length of the top blast hit. 
Chimeric transcripts are going to be a big problem for this approach, but I can check for 
them independently (by comparing cds and exon length) and break them up if necessary before 
doing the expression analysis

I am going to modify this to parse start, end and strand info from the fasta header of aa 
sequences from OrfPredictor (http://proteomics.ysu.edu/tools/OrfPredictor.html)
"""

import sys  
import getopt 
import csv
from os import path
import cPickle as pickle
from Heathpy import flatten_GTF
from Bio import SeqIO

def main(argv):
  blastfilename = ''
  gtf_filename = ''
  aa_filename = ''
  try:
    opts, args = getopt.getopt(argv,"hb:g:a:",["blast=","gtf=","aa"])
  except getopt.GetoptError:
    print 'Type Blast2GTF.py -h for options'
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print 'Blast2GTF.py -b <blastfile> -g <GTF_file> -a <AA_file>'
       sys.exit()
    elif opt in ("-g", "--gtf"):
       gtf_filename = arg
    elif opt in ("-b", "--blast"):
       blastfilename = arg
    elif opt in ("-a", "--aa"):
       aa_filename = arg

  previous_query = 'init'
  cds = {
  'source': '1kp',
  'feature': 'CDS',
  'frame': 1
  }

  #make a list of names to ensure that there are no duplicates
  name_list = {}
  name_num = 1
  
  #read dictionary of cluster membership
  cluster_info = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "RefSeq", "SeqClusters.p")
  seq_groups = pickle.load( open( cluster_info, "rb" ) )
  
  #Read through AA sequences and parse header info to get start, end and strand values
  orf_predictions = {}
  if aa_filename:
    for seq_record in SeqIO.parse(aa_filename, "fasta"):
      orf_predictions[seq_record.id] = Orf_stats(seq_record.description)    
  
  with open(gtf_filename, 'wb') as outfile:
    gtf_writer = csv.writer(outfile, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
    with open(blastfilename, 'rU') as infile:
      reader=csv.reader(infile,delimiter='\t')
      for row in reader:
        qseqid, qlen, sacc, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore = row
        qlen = int(qlen)
        slen = int(slen)
        pident = float(pident)
        length = int(length)
        mismatch = int(mismatch)
        gapopen = int(gapopen)
        qstart = int(qstart)
        qend = int(qend)
        qframe =int(qframe)
        sstart = int(sstart)
        send = int(send)
        sframe = int(sframe)
        evalue = float(evalue)
        bitscore = float(bitscore)
                 
        #skip secondary hits from same query
        if qseqid == previous_query: #skip secondary hits
          continue 

        previous_query = qseqid
    
        #Get info for exon feature
        cds['seqname'] = qseqid
        cds['score'] = bitscore
        if qseqid in orf_predictions:                          #Use data from OrfPredictor if present
          cds['start'] = orf_predictions[qseqid].start
          cds['end'] = orf_predictions[qseqid].end
          cds['strand'] = orf_predictions[qseqid].strand
        elif sstart < send and qstart < qend:                 #I'm assuming there's no reason why both would be reversed
          cds['strand'] = '+'
          cds['start'] = qstart
          cds['end'] = qend
        else:
          cds['strand'] = '-'
          cds['start'] = qend
          cds['end'] = qstart   
        if sacc in seq_groups:
          name = "%s_c%s_0" % (qseqid[0:qseqid.find('comp')], seq_groups[sacc].split("_")[1])
        else:
          name =  sacc + "_0"                       
        while name in name_list:
          name = name.split("_")[0:-1] + [str(name_num)]
          name = "_".join(name)
          name_num = name_num + 1
        name_list[name] = 1
        name_num = 1
        cds['gene_id'] = name
        cds['transcript_id'] =  cds['gene_id'] + '.1'
        gtf_writer.writerow(flatten_GTF(cds))

class Orf_stats:
  def __init__(self, header):
    attributes = header.split("\t")
    self.start = attributes[2]
    self.end = attributes[3]
    if '+' in attributes[1]:
      self.strand = '+'
    else:
      self.strand = '-'
      
if __name__ == "__main__":
   main(sys.argv[1:])


