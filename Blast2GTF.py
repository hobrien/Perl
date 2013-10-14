#!/usr/local/bin/python

"""This is being rewritten to make a single cds feature the length of the top blast hit. 
Chimeric transcripts are going to be a big problem for this approach, but I can check for 
them independently (by comparing cds and exon length) and break them up if necessary before 
doing the expression analysis

I am going to modify this to parse start, end and strand info from the fasta header of aa 
sequences from OrfPredictor (http://proteomics.ysu.edu/tools/OrfPredictor.html)

The new plan (10 Oct 2013) is to record top non-overlapping blast hits as "blast_hit" features
and the surrounding ORFs as CDS features. Exon info can be merged from cufflinks later.

If there are multiple hits to the same target sequence, a warning will be issued.

It there are multiple hits to different targets and they don't overlap, all will be written
"""

import sys, warnings , logging 
import getopt 
import csv
from os import path
import cPickle as pickle
from Heathpy import flatten_GTF, get_orf_coords, warning_on_one_line
from Bio import SeqIO

def main(argv):
  blastfilename = ''
  gtf_filename = ''
  seq_file = ''
  try:
    opts, args = getopt.getopt(argv,"hb:g:s:",["blast=","gtf=","sequence="])
  except getopt.GetoptError:
    print 'Type Blast2GTF.py -h for options'
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print 'Blast2GTF.py -b <blastfile> -g <GTF_file> -s <seqfile>'
       sys.exit()
    elif opt in ("-g", "--gtf"):
       gtf_filename = arg
    elif opt in ("-b", "--blast"):
       blastfilename = arg
    elif opt in ("-s", "--sequence"):
       seq_file = arg
  
  LOG_FILENAME = '.'.join(gtf_filename.split('.')[:-1] + ['log'])
  logging.basicConfig(filename=LOG_FILENAME,level=logging.WARNING)

  previous_query = 'init'
  previous_hits = {}
  frameshift_warn = 0
    
  #make a list of names to ensure that there are no duplicates
  name_list = {}
  name_num = 1
  
  #read dictionary of cluster membership
  cluster_info = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "RefSeq", "SeqClusters.p")
  seq_groups = pickle.load( open( cluster_info, "rb" ) )
  seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
  
  with open(gtf_filename, 'wb') as outfile:
    gtf_writer = csv.writer(outfile, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
    with open(blastfilename, 'rU') as infile:
      reader=csv.reader(infile,delimiter='\t')
      for row in reader:
        qseqid, qlen, sacc, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore = row
        new_hit = {
          'qseqid': qseqid,
          'qlen': int(qlen),
          'sacc': sacc,
          'slen': int(slen),
          'pident': float(pident),
          'length': int(length),
          'mismatch': int(mismatch),
          'gapopen': int(gapopen),
          'qstart': int(qstart),
          'qend': int(qend),
          'qframe': int(qframe),
          'sstart': int(sstart),
          'send': int(send),
          'sframe': int(sframe),
          'evalue': float(evalue),
          'bitscore': float(bitscore)
        }
        #Add strand information and reverse coordinates if on negative strand
        new_hit['strand'] = '+'
        if new_hit['qend'] < new_hit['qstart']:
          new_hit['strand'] = '-'
          new_hit['qstart'] = int(qend)
          new_hit['qend'] = int(qstart)
        if new_hit['send'] < new_hit['sstart']:
          new_hit['strand'] = '-'
          new_hit['sstart'] = int(send)
          new_hit['send'] = int(sstart)
                  
        #Make list of all non-overlapping hits, printing a warning if there are multiple hits to the same sequence  
        if new_hit['qseqid'] == previous_query:
          if new_hit['sacc'] in previous_hits:
            if frameshift_warn == 0:
              logging.warning("%s has multiple hits to %s" % (new_hit['qseqid'], new_hit['sacc']))
              frameshift_warn = 1
          else:
            overlap = 0
            for hit in previous_hits.keys():
              overlap = max(overlap, hit_overlap(previous_hits[hit], new_hit))
            if not overlap:
              previous_hits[new_hit['sacc']] = new_hit

        else:
          for hit in previous_hits.values():
            #write info about blast hit to GTF
            feature = {'source': '1kp',
                       'feature': 'blast_hit',
                       'frame': 1,
                       'seqname': hit['qseqid'],
                       'score': hit['bitscore'],
                       'strand': hit['strand'],
                       'start': hit['qstart'],
                       'end': hit['qend']   
                      }

            if hit['sacc'] in seq_groups:
              name = "%s_c%s_0" % (qseqid[0:qseqid.find('comp')], seq_groups[hit['sacc']].split("_")[1])
            else:
              name =  hit['sacc'] + "_0"                       
            while name in name_list:
              name = name.split("_")[0:-1] + [str(name_num)]
              name = "_".join(name)
              name_num = name_num + 1
            name_list[name] = 1
            name_num = 1
            feature['gene_id'] = name
            feature['transcript_id'] =  feature['gene_id'] + '.1'
            gtf_writer.writerow(flatten_GTF(feature))
            
            #write info about ORF containing blast hit to file
            feature['feature'] = 'CDS'
            feature['score'] = '.'
            seq = seq_dict[feature['seqname']]
            if feature['strand'] == '+':
              (feature['start'], feature['end']) = get_orf_coords(seq, feature['start'], feature['end'])
            else:
              (feature['start'], feature['end']) = get_orf_coords(seq.reverse_complement(), len(seq) - feature['end'] + 1, len(seq) - feature['start'] + 1)              
            gtf_writer.writerow(flatten_GTF(feature))
            
          previous_hits = {new_hit['sacc']: new_hit }
          previous_query = qseqid
          frameshift_warn = 0

def hit_overlap(hit1, hit2):
  overlap = 0
  if hit1['qstart'] > hit2['qstart'] and hit1['qstart'] < hit2['qend']: #hit1 starts within hit2
    overlap = 1
  if hit1['qend'] > hit2['qstart'] and hit1['qend'] < hit2['qend']:     #hit1 ends within hit2
    overlap = 1
  if hit1['qstart'] < hit2['qstart'] and hit1['qend'] > hit2['qend']:  #hit1 starts before hit2 and ends after
    overlap = 1
  return overlap  
  
      
if __name__ == "__main__":
   main(sys.argv[1:])


