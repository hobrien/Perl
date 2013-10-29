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

Now (22 Oct) I'm going to refactor it to work with MySQL (all of the data is now in a MySQL db
and it's actually indexed so the lookups are fast.

This is also going to involve a bit more of a rewrite because the sequences are in the db

"""

import sys, warnings , logging 
import getopt 
import csv
import MySQLdb as mdb
from os import path
import cPickle as pickle
from Heathpy import flatten_GTF, get_orf_coords, warning_on_one_line, make_db
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(argv):
  gtf_filename = ''
  species = ''
  seqfilename = ''
  usage = 'Blast2GTF.py -n <species_name> -g <gtf_outfile> -s <seq_file>'
  try:
    opts, args = getopt.getopt(argv,"hn:g:s:",["name=", "gtf=", "seqfile="])
    if not opts:
      raise getopt.GetoptError('no opts')
  except getopt.GetoptError:
    print usage
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print usage
       sys.exit()
    elif opt in ("-n", "--name"):
       species = arg
    elif opt in ("-s", "--seqfile"):
       seqfilename = arg
    elif opt in ("-g", "--gtf"):
       gtf_filename = arg
  
  
  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
    #make a list of names to ensure that there are no duplicates
    name_list = {}
      
    #read dictionary of cluster membership (It would be far better to put this info into the db, but this is working, so I won't mess with it
    cluster_info = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "RefSeq", "SeqClusters.p")
    seq_groups = pickle.load( open( cluster_info, "rb" ) )
    
    if gtf_filename:
      outfile = open(gtf_filename, 'wb')
      gtf_writer = csv.writer(outfile, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
    
    if seqfilename:
      seqfile = open(seqfilename, 'wb')
    cur = con.cursor()
    cur.execute("SELECT a.seqid, b.sequence FROM Species a, Ortholog_groups b WHERE a.seqid = b.seqid AND a.species= %s", (species))
    #cur.execute("SELECT seqid, sequence FROM Ortholog_groups b WHERE seqid = 'UNCcomp100025_c0_seq1'")
    rows = cur.fetchall()
    for (seqid, seq) in rows:
      try:
        seq_record = SeqRecord(Seq(seq))
      except TypeError:
        warnings.warn("Can't create seq object from %s for %s" % (seq, seqid))
        continue
      if gtf_filename:       #Write blast info to GTF file
        hit_list = {}
        frameshift_warn = 0
        cur.execute("SELECT * FROM BLAST WHERE qseqid=%s", (seqid,))
        for (id, qseqid, qlen, sacc, slen, pident, length, mismatch, gapopen, qstart, qend, qframe, sstart, send, sframe, evalue, bitscore, strand, hitnum) in cur.fetchall():
          feature = {
            'source': '1kp',
            'feature': 'blast_hit',
            'frame': '.',
            'seqname': qseqid,
            'score': float(bitscore),
            'start': int(qstart),
            'end': int(qend)
                   }
          #Add strand information and reverse coordinates if on negative strand
          if strand == '1':
            feature['strand'] = '+'
          elif strand == '0':
            feature['strand'] = '-'
          else:
            sys.exit("Strand %s not recognized" % strand)        
          #Make list of all non-overlapping hits, printing a warning if there are multiple hits to the same sequence
          if sacc in hit_list:
            if frameshift_warn == 0:
              warnings.warn("%s has multiple hits to %s" % (qseqid, sacc))
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
              name_num = 1                         
              while name in name_list:
                name = name.split("_")[0:-1] + [str(name_num)]
                name = "_".join(name)
                name_num = name_num + 1
              name_list[name] = 1
              feature['gene_id'] = name
              feature['transcript_id'] =  feature['gene_id'] + '.1'
              hit_list[sacc] = feature

        for feature in hit_list.values():
          #write info about blast hit to GTF
          #gtf_writer.writerow(flatten_GTF(feature))
            
          #write info about ORF containing blast hit to file
          feature['feature'] = 'CDS'
          feature['score'] = '.'
          if feature['strand'] == '+':
            (feature['start'], feature['end']) = get_orf_coords(seq_record, feature['start'], feature['end'])
            feature['frame'] = feature['start'] % 3 + 1
            #print "feature start, feature end = %s, %s" % ( feature['start'], feature['end'] )
            note = orf_integrity(Seq(seq[feature['start']-1:feature['end']]))
            if note: feature['note'] = note
          else:
            (orf_start, orf_end) = get_orf_coords(seq_record.reverse_complement(), len(seq_record) - qend + 1, len(seq_record) - qstart + 1)    
            (feature['start'], feature['end']) = (len(seq_record) - orf_end + 1, len(seq_record) - orf_start + 1)            
            feature['frame'] = ( len(seq_record) - feature['end'] + 1 ) % 3 + 1
            orf_rev = Seq(seq[feature['start']-1:feature['end']] )          
            note = orf_integrity(orf_rev.reverse_complement())
            if note: feature['note'] = note
          gtf_writer.writerow(flatten_GTF(feature))

      if seqfilename:
        seqfile.write(">%s\n%s\n" % (seqid, seq))

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

def orf_integrity(seq):
  orf = ''
  aa_seq = seq.translate()
  if aa_seq[0] != 'M':
    orf = "C-terminal fragment"
  end = aa_seq[-1]
  if end != '*':
    if orf:
      orf = "Internal fragment"
    else:
      orf = "N-terminal fragment"
  if seq[:-1].find("*") > -1:
    orf = "Pseudogene"
  return orf

if __name__ == "__main__":
   main(sys.argv[1:])


