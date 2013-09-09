#!/opt/local/bin/python


import sys
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

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



def smart_consensus(aln, threshold = .3, ambiguous = "N", gap = "-",
                  require_multiple = 0, consensus_alpha = "IUPACUnambiguousDNA"): 
    """stolen from biopython but modified to ignore ambiguities (while including gaps)
       this is useful in cases where sequences vary in length rather than in actual bases
    
       gap can be set to "" to produce a continuous sequence which will not match the original alignment
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

        if require_multiple and num_atoms == 1: 
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
  
if __name__ == "__main__":
   aln = AlignIO.read(sys.argv[1], sys.argv[2])  #call with 2 arguments: filename and format
   print smart_consensus(pad_ends(aln))
