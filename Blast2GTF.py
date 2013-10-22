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

Now (15 Oct) I'm going to refactor this to work with a SQLite database of blast results.

"""

import sys, warnings , logging 
import getopt 
import csv
import sqlite3
from os import path
import cPickle as pickle
from Heathpy import flatten_GTF, get_orf_coords, warning_on_one_line, make_db
from Bio import SeqIO

def main(argv):
  blastfilename = ''
  gtf_filename = ''
  seqfilename = ''
  usage = 'Blast2GTF.py -b <blastfile> -g <GTF_file> -s <seqfile>'
  try:
    opts, args = getopt.getopt(argv,"hb:g:s:",["blast=","gtf=","sequence="])
    if not opts:
      raise getopt.GetoptError('no opts')
  except getopt.GetoptError:
    print usage
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print usage
       sys.exit()
    elif opt in ("-g", "--gtf"):
       gtf_filename = arg
    elif opt in ("-b", "--blast"):
       blastfilename = arg
    elif opt in ("-s", "--sequence"):
       seqfilename = arg
  
  LOG_FILENAME = '.'.join(gtf_filename.split('.')[:-1] + ['log'])
  logging.basicConfig(filename=LOG_FILENAME,level=logging.WARNING)

  #Build database of blast results (if one doesn't already exist)
  db_name = blastfilename + '.db'
  if not path.exists(db_name):
    make_db(blastfilename, db_name)

  
    
  #make a list of names to ensure that there are no duplicates
  name_list = {}
  name_num = 1
  
  #read dictionary of cluster membership
  cluster_info = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "RefSeq", "SeqClusters.p")
  seq_groups = pickle.load( open( cluster_info, "rb" ) )

  conn = sqlite3.connect(db_name)
  c = conn.cursor()

  outfile = open(gtf_filename, 'wb')
  gtf_writer = csv.writer(outfile, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
  
  for seq_record in SeqIO.parse(seqfilename, "fasta"):
    hit_list = {}
    frameshift_warn = 0
    for row in c.execute('SELECT * FROM hits WHERE qseqid=?', (seq_record.id,)): 
      qseqid, qlen, sacc, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore, strand, hitnum = row
      feature = {
        'source': '1kp',
        'feature': 'blast_hit',
        'frame': 1,
        'seqname': qseqid,
        'score': float(bitscore),
        'start': int(qstart),
        'end': int(qend)
      }
      #Add strand information and reverse coordinates if on negative strand
      if strand:
        feature['strand'] = '+'
      else:
        feature['strand'] = '-'
                  
      #Make list of all non-overlapping hits, printing a warning if there are multiple hits to the same sequence  
      if sacc in hit_list:
        if frameshift_warn == 0:
          logging.warning("%s has multiple hits to %s" % (qseqid, sacc))
          frameshift_warn = 1
      else:
        overlap = 0
        for hit in hit_list.keys():
          overlap = max(overlap, hit_overlap(hit_list[hit], feature))
        if not overlap:
          if sacc in seq_groups:
            name = "%s_c%s_0" % (qseqid[0:qseqid.find('comp')], seq_groups[sacc].split("_")[1])
          else:
            name = sacc + "_0"                       
          while name in name_list:
            name = name.split("_")[0:-1] + [str(name_num)]
            name = "_".join(name)
            name_num = name_num + 1
          name_list[name] = 1
          name_num = 1
          feature['gene_id'] = name
          feature['transcript_id'] =  feature['gene_id'] + '.1'
          hit_list[sacc] = feature

    for feature in hit_list.values():
      #write info about blast hit to GTF
      gtf_writer.writerow(flatten_GTF(feature))
            
      #write info about ORF containing blast hit to file
      feature['feature'] = 'CDS'
      feature['score'] = '.'
      if feature['strand'] == '+':
        (feature['start'], feature['end']) = get_orf_coords(seq_record, feature['start'], feature['end'])
      else:
        (feature['start'], feature['end']) = get_orf_coords(seq_record.reverse_complement(), len(seq_record) - feature['end'] + 1, len(seq_record) - feature['start'] + 1)              
      gtf_writer.writerow(flatten_GTF(feature))


def hit_overlap(hit1, hit2):
  overlap = 0
  #print "Comparing hit1 (%i - %i) to hit2 (%i - %i)" % (hit1['start'], hit1['end'], hit2['start'], hit2['end'])
  if hit1['start'] >= hit2['start'] and hit1['start'] <= hit2['end']: #hit1 starts within hit2
    overlap = 1
  if hit1['end'] >= hit2['start'] and hit1['end'] <= hit2['end']:     #hit1 ends within hit2
    overlap = 1
  if hit1['start'] <= hit2['start'] and hit1['end'] >= hit2['end']:  #hit1 starts before hit2 and ends after
    overlap = 1
  return overlap  
  
      
if __name__ == "__main__":
   main(sys.argv[1:])


