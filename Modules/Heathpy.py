#!/opt/local/bin/python


import sys, warnings, string
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

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
    fields.append(feature['score'])
    del feature['score']
  except KeyError:
    sys.exit("no score tag")
  try:
    fields.append(feature['strand'])
    del feature['strand']
  except KeyError:
    sys.exit("no strand tag")
  try:
    fields.append(feature['frame'])
    del feature['frame']
  except KeyError:
    sys.exit("no frame tag")
  attributes = ""
  keys = feature.keys().sort()
  for key in sorted(feature, key=feature.get):
    attributes = attributes + '%s "%s"; ' % (key, feature[key])
  fields.append(attributes)
  return fields

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
  CDS = seq[cds_start:cds_end+1]

  return CDS 
   
def get_orf_coords(seq, blast_start, blast_end):
  warnings.formatwarning = warning_on_one_line
  #determine the number of bp that must be trimmed from the start of the sequence to give correct frame
  frame = blast_start % 3 - 1
  if frame == -1:
    frame = 2

  #convert DNA coordinates to AA translation coordinates
  tr_start = ( blast_start - frame ) / 3
  tr_end = ( blast_end - frame -2 ) / 3  #subtract 2 extra bp because blast reports the last bp of the codon. We want the first

  #Trim sequence to correct frame and remove partial codons (up to 2 bp from each end)
  to_translate = seq[frame:]          #trim start
  while len(to_translate) % 3:        #trim partial codons from end
    to_translate = to_translate[:-1]

  aa_tr = to_translate.seq.translate()  #Translate sequence

  #Check for stop codons within homologous region (nonsense mutation)
  if string.find(aa_tr, "*", tr_start, tr_end - 1) != -1:
    warnings.warn("%s: Sequence contains a stop codon within the homologous region (likely due to a nonsense mutation)" % seq.id)

  #Determine the postion of the first in-frame stop after (or at the end of) the homologous region  
  aa_end = string.find(aa_tr, "*", tr_end - 1)  #This includes the last AA of the blast hit

  if aa_end == -1:                              #No stop codon found. ORF extends beyond the contig
    warnings.warn("%s: No stop codon. It is likely that the ORF extends beyond the contig (N-terminal fragement)" % seq.id)
    aa_end = len(aa_tr)

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
  cds_start = aa_start*3+frame
  cds_end = aa_end*3+frame+2      #Add extra two bases to include all of final codon
  return (cds_start, cds_end)
 
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