#!/opt/local/bin/python


import sys, warnings, string
import sqlite3
import csv
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

def make_db(blastfile, db_name):
  conn = sqlite3.connect(db_name)
  c = conn.cursor()
  #strand: 1 = pos, 0 = neg
  c.execute('''CREATE TABLE hits
                (qseqid text, qlen int, sacc text, slen int, pident read, length int, mismatch int, gapopen int, qstart int, qend int, qframe int, sstart int, send int, sframe int, evalue float, bitscore float, strand bit, hitnum int)''')
  with open(blastfile, 'rU') as infile:
    reader=csv.reader(infile,delimiter='\t')
    prev = 'init'
    hit_num = 0
    for row in reader:
      row.append(1)  #add strand info to row
      if int(row[8]) > int(row[9]):               #qstart > qend
        (row[8], row[9]) = (row[9], row[8])       #reverse qstart and qend
        row[-1] = 0                               #change strand to neg
      if int(row[11]) > int(row[12]):             #sstart > send
        (row[11], row[12]) = (row[12], row[11])   #reverse sstart and send
        row[-1] = 0                               #change strand to neg
      if row[0] == prev:
        hit_num += 1
      else:
        prev = row[0]
        hit_num = 1
      row.append(hit_num)                                      #add count num to row
      c.execute("INSERT INTO hits VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,?)", row)
    conn.commit()
  conn.close()

def parse_GTF (fields):
  tags = {}
  for attribute in fields[8].split(";")[:-1]:
    attribute = attribute.strip()
    tags[attribute.split(" ")[0]] = " ".join(attribute.split(" ")[1:]).replace('"','')
  try:
    tags['frame'] = int(fields[7])
  except ValueError:
    if fields[7] == '.':
      tags['frame'] = 'NULL'
    else:
      sys.exit("frame %s not recognized. Must be 1, 2, 3 or ." % fields[5]) 
  if fields[6] == '-' or fields[6] == 0 or fields[6] == -1:
    tags['strand'] = 0
  elif fields[6] == '+' or fields[6] == 1:
    tags['strand'] = 1
  elif fields[6] == '.':
    tags['strand'] = 'NULL'
  else:  
    sys.exit("strand %s not recognized. Must one of +, -, 1, -1 or 0" % fields[6])      
  try:
    tags['score'] = float(fields[5])
  except ValueError:
    if fields[5] == '.':
      tags['score'] = 'NULL'
    else:
      sys.exit("score %s not recognized. Must be a number" % fields[5])
  try:
    tags['end'] = int(fields[4])
  except ValueError:
    sys.exit("score %s not recognized. Must be a positive integer" % fields[4])
  try:
    tags['start'] = int(fields[3])
  except ValueError:
    sys.exit("score %s not recognized. Must be a positive integer" % fields[3])
  tags['feature'] = fields[2]
  tags['source'] = fields[1]
  tags['seqid'] = fields[0]
  return tags



def flatten_GTF(input):
  feature = input.copy()
  fields = []
  try:
    fields.append(feature['seqname'])
    del feature['seqname']
  except KeyError:
    sys.exit("no seqname tag")
  try:
    fields.append(feature['source'])
    del feature['source']
  except KeyError:
    sys.exit("no source tag")
  try:
    fields.append(feature['feature'])
    del feature['feature']
  except KeyError:
    sys.exit("no feature tag")
  try:
    fields.append(feature['start'])
    del feature['start']
  except KeyError:
    sys.exit("no start tag")
  try:
    fields.append(feature['end'])
    del feature['end']
  except KeyError:
    sys.exit("no end tag")
  try:
    if feature['score']:
      fields.append(feature['score'])
    else:
      fields.append('.')
    
    del feature['score']
  except KeyError:
    sys.exit("no score tag")
  try:
    if int(feature['strand']) == 0:
      fields.append('-')
    elif int(feature['strand']) == 1:
      fields.append('+')
    elif feature['strand'] == '.':
      fields.append('.')
    else:
      sys.exit("strand %s not recognized" % feature['strand'])
    del feature['strand']
  except KeyError:
    sys.exit("no strand tag")
  try:
    fields.append(feature['frame'])
    del feature['frame']
  except KeyError:
    sys.exit("no frame tag")
  attributes = ""
  keys = feature.keys()
  keys.sort()
  #for key in sorted(feature, key=feature.get):
  for key in keys:
    attributes = attributes + '%s "%s"; ' % (key, feature[key])
  fields.append(attributes)
  return "\t".join(map(str, fields))

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    import warnings
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)

def smart_consensus(aln, threshold = .3, ambiguous = "N", gap = "-",
                  require_multiple = 0, consensus_alpha = "IUPACUnambiguousDNA"): 
    """stolen from biopython but modified to ignore ambiguities (while including gaps)
       this is useful in cases where sequences vary in length rather than in actual bases
    
       gap can be set to "" to produce a continuous sequence which will not match the original alignment
       
       this is failing in cases where all sequences have gaps at a position. 
    """
    consensus = '' 

    # find the length of the consensus we are creating 
    con_len = aln.get_alignment_length() 

    # go through each seq item 
    for n in range(con_len): 
        # keep track of the counts of the different atoms we get 
        atom_dict = {} 
        num_atoms = 0 

        for record in aln._records: 
            # make sure we haven't run past the end of any sequences 
            # if they are of different lengths 
            if n < len(record.seq):
                if record.seq[n] != ambiguous:  #skip sequences that have ambiguity
                    if record.seq[n] not in atom_dict: 
                        atom_dict[record.seq[n]] = 1 
                    else: 
                        atom_dict[record.seq[n]] += 1 

                    num_atoms += 1 

        max_atoms = [] 
        max_size = 0 

        for atom in atom_dict: 
            if atom_dict[atom] > max_size: 
                max_atoms = [atom] 
                max_size = atom_dict[atom] 
            elif atom_dict[atom] == max_size:  #if gap character and a base are equally common, use gap
                if atom == gap:
                    max_atoms = [atom]
                elif gap not in max_atoms:
                    max_atoms.append(atom) 
        if num_atoms == 0:                    #if all sequences have an N at a position (indicates that all sequences start or stop before/after this position)
            consensus += gap
        elif require_multiple and num_atoms == 1: 
            consensus += ambiguous 
        elif (len(max_atoms) == 1) and ((float(max_size)/float(num_atoms)) 
                                     >= threshold): 
            consensus += max_atoms[0] 
        else: 
            consensus += ambiguous 


    return Seq(consensus, consensus_alpha) 


def pad_ends(aln, ambig = "N", gap = "-"):
  """Replaces gaps at the start or end of aligned sequences with ambiguities"""
  padded_seqs = []
  for seq in aln:
    #print seq.seq
    start = 0
    while start < len(seq) and seq[start] == gap:
      start = start + 1
      
    end = len(seq) - 1
    while end > 0 and seq[end] == gap:
      end =end - 1
    end = end + 1
    end_pad = len(seq) - end
    #print "seq length: %s, start: %s, end: %s, end_pad: %s" % (len(seq), start, end, end_pad)
    seq = ambig*start + seq[start:end] + ambig*(end_pad)
    #print seq.seq
    padded_seqs.append(seq)
  return padded_seqs

def max_key (dict):
  max_keys = [] 
  max_value = 0 
  for key in dict: 
    if dict[key] > max_value: 
      max_keys = [key] 
      max_value = dict[key] 
    elif dict[key] == max_value: 
      max_keys.append(key) 
  return max_keys

def get_orf(seq, blast_start, blast_end):
  (cds_start, cds_end) = get_orf_coords(seq, blast_start, blast_end)
  CDS = seq[cds_start - 1:cds_end]  #Convert from 1 based counting

  return CDS 
   
def get_orf_coords(seq, blast_start, blast_end):
  warnings.formatwarning = warning_on_one_line
  #determine the number of bp that must be trimmed from the start of the sequence to give correct frame
  start_frame = blast_start % 3 - 1
  if start_frame == -1:
    start_frame = 2


  #Trim sequence to correct frame and remove partial codons (up to 2 bp from each end)
  to_translate = seq[start_frame:]          #trim start
  while len(to_translate) % 3:        #trim partial codons from end
    to_translate = to_translate[:-1]

  aa_tr = to_translate.seq.translate()  #Translate sequence

  #If there is a frameshift, I need to modify the sequence then retranslate so that the end is in frame
  end_frame = (blast_end -2) % 3 - 1 #subtract 2 extra bp because blast reports the last bp of the codon. We want the first
  if end_frame == -1:
    end_frame = 2
  if end_frame == start_frame:
    aa_tr_end = aa_tr
  else:
    to_end = seq[end_frame:]
    while len(to_end) % 3:        #trim partial codons from end
      to_end = to_end[:-1]
    aa_tr_end = to_end.seq.translate()  #Translate sequence
    
  #convert DNA coordinates to AA translation coordinates
  tr_start = ( blast_start - start_frame ) / 3
  tr_end = ( blast_end - end_frame -2 ) / 3  #subtract 2 extra bp because blast reports the last bp of the codon. We want the first

  #Check for stop codons within homologous region (nonsense mutation)
  if string.find(aa_tr, "*", tr_start, tr_end - 1) != -1:
    warnings.warn("%s: Sequence contains a stop codon within the homologous region (likely due to a nonsense mutation)" % seq.id)

  #Determine the postion of the first in-frame stop after (or at the end of) the homologous region
  
  aa_end = string.find(aa_tr_end, "*", tr_end - 1)  #This includes the last AA of the blast hit
  if aa_end == -1:                              #No stop codon found. ORF extends beyond the contig
    warnings.warn("%s: No stop codon. It is likely that the ORF extends beyond the contig (N-terminal fragement)" % seq.id)
    aa_end = len(aa_tr) - 1
    

  #determine position of LAST in-frame stop before homologous region (start codon must occur after this position)
  minimum_start = string.rfind(aa_tr, "*", 0, tr_start)  
  if minimum_start == -1:          #no stop codon before homologous region
    minimum_start = 0

  #find first met after last in-frame start before (or at) start of homologous region
  # (could also look for last met if we want to minimize the size of non-homologous coding sequence)
  aa_start = string.find(aa_tr, "M", minimum_start, tr_start + 1) #include first AA of blast hit

  if aa_start == -1:              #No start codon found. Either start ORF extends beyond
    if minimum_start == 0:
      warnings.warn("%s: No start codon.  It is likely that the ORF extends beyond the contig (C-terminal fragement).\n\t\t\t\t\t Setting ORF start to the start of the homologous region" % seq.id)
    else:
      warnings.warn("%s: Stop codon between start codon and homologous region (likely due to a nonsense mutation in non-homologous region).\n\t\t\t\t\t Setting ORF start to the start of the homologous region" % seq.id)
    aa_start = tr_start

  #Convert ORF coordinates back to DNA coordinates ( adding frame so numbers correspond to original seq.
  #Could also use trimmed DNA seq, in which case this is unnecessary
  cds_start = aa_start*3+start_frame
  cds_end = (aa_end*3)+end_frame+2       #Add extra two bases to include all of final codon
  return (cds_start+1, cds_end+1)  #Convert to one-based counting
 
if __name__ == "__main__":
  #aln = AlignIO.read(sys.argv[1], sys.argv[2])  #call with 2 arguments: filename and format
  #print smart_consensus(pad_ends(aln))
  seq = SeqRecord(Seq("TTTAGTTTTTTATGTTTTTTTTTTAGTTTTAG", IUPAC.unambiguous_dna), id='test_seq')
  print "Testing get_orf with intact ORF:"
  print  get_orf(seq, 15, 23).seq
  print "Testing get_orf with pseudogene:"
  print get_orf(seq, 15, 32).seq
  print "Testing get_orf with no start codon:"
  print get_orf(seq, 6, 23).seq  
  print "Testing get_orf with no stop codon:"
  print get_orf(seq[:23], 15, 21).seq  
  print "Testing get_orf with stop codon between start and homologous region:"
  print get_orf(seq, 27, 32).seq
  print "Testting get_orf with frameshift:"
  seq = SeqRecord(Seq("TTTAGTTTTTTATGTTTTTTTTTTTAGTTTTAG", IUPAC.unambiguous_dna), id='test_seq')
  print  get_orf(seq, 15, 21).seq
  