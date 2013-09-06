#!/opt/local/bin/python
"""

MakeConsesus.py 5 Sept 2013

SYNOPSIS

MakeConsensus.py -t <treefile> -a <alignmentfile> -c <consensusfile> -l <logfile>

Options:
 -t treefile
 -a alignment file (relaxed Phylip format)
 -c outfile for consensus sequences (appended)
 -l logfile to record summary of which sequences are represented by the consensus (appended)

DESCRIPTION
Takes a alignment file and resultant treefile and looks for monophyletic groups of sequences
from a single species, then makes consensus sequences from these groups

Consensus sequence is appended to the specified file and a summary of which sequences it
represents are appended to the specified logfile

NOTES
Consensus considers INTERNAL gaps as a valid character state (ie; if a plurality of sequences
have a gap at a position, it is considered a gap in the consensus), but ignores gaps
before/after the sequence

Consensus sequences are named after the most common comp number among the included sequences
(in case of a tie, one is picked arbitrarily). This means that there is a possibility that
if members of the same comp are found in multiple clusters, there may end up being multiple
consensus sequences with the same name. This should be uncommon, but is something to watch
out for.

At the moment, this will only make a consensus if ALL of the sequences fro a species are
monophyletic. This will have to be improved in the future.

AUTHOR

 Heath E. O'Brien (heath.obrien-at-gmail-dot-com)
"""

import sys  
import getopt 
from os import path, system
from ete2 import Tree  
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Heathpy import smart_consensus, max_key, pad_ends

def main(argv):
   treefile = ''
   alnmentfile = ''
   consfilename = ''
   logfilename = ''
   try:
      opts, args = getopt.getopt(argv,"ht:a:c:l:",["tree=","alignment=","consensus=","log="])
   except getopt.GetoptError:
      print 'MakeConsensus.py -t <treefile> -a <alignmentfile> -c <consensusfile> -l <logfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'MakeConsensus.py -t <treefile> -a <alignmentfile> -c <consensusfile> -l <logfile>'
         sys.exit()
      elif opt in ("-t", "--tree"):
         treefile = arg
      elif opt in ("-a", "--alignment"):
         alnfile = arg
      elif opt in ("-c", "--consensus"):
         consfilename = arg
      elif opt in ("-l", "--log"):
         logfilename = arg
   
   
   #read in tree and root by midpoint
   t = Tree(treefile)
   R = t.get_midpoint_outgroup()
   t.set_outgroup(R)
   
   #read in alignment and create index of sequence ids to access individual sequences by id
   aln = AlignIO.read(alnfile, "phylip-relaxed")
   aln_index = {}
   seq_num = 0
   for seq in aln:
     aln_index[seq.id] = seq_num
     seq_num = seq_num +1
   
   #loop through list of species andprint consensus if seqs are monophyletic
   species_list = ['UNC', 'MOEL', 'WILD', 'KRAUS'] #, 'REF']
   for species in species_list:
     leaves = GetLeaves(t, species)
     if len(leaves) > 1 and MonoTest(t, species) == 0:
       
       #make list of sequences from species
       species_seqs = []
       for leaf in leaves:
         if 'kraussiana' not in leaf.name:  #exclude reference sequence from consensus
           species_seqs.append(aln[aln_index[leaf.name]])

       #create consensus sequence from sequence list
       species_aln = MultipleSeqAlignment(species_seqs)
       consensus = smart_consensus(MultipleSeqAlignment(pad_ends(species_aln))) #use pad_ends to convert end gaps to ambiguous
       consensus = pad_ends([consensus], '-', 'N')[0]  #convert back to gap characters (need to convert seq to a list)
       name = pick_name(leaves) + "_cons"
       consensus =  SeqRecord(Seq(str(consensus).replace('-','')),
                              id=name, 
                              description=name)
                              
       #write consensus sequence to consensus file
       consfile = open(consfilename, "a")
       consfile.write(consensus.format("fasta"))
       consfile.close()
       
       #Remove individual sequences from the cluster and add consensus
       seqfilename = alnfile.split(".")[0] + ".fa"
       system("RemoveSeqs.py -s %s -g %s" % (seqfilename, species))
       seqfile = open(seqfilename, "a")
       seqfile.write(consensus.format("fasta"))
       seqfile.close()
       
       #write summary of seqs represented by cons to log file
       cluster_num = path.split(alnfile)[1].split(".")[0] #get cluster number
       logfile = open(logfilename, "a")
       for seq in species_seqs:
         logfile.write("%s\n" % ", ".join([seq.id, name, cluster_num]))
       logfile.close()
 
def pick_name(leaves):
  names = {}
  for leaf in leaves:
    name = leaf.name.split("_")[0]
    if name in names:
      names[name] = names[name] + 1
    else:
      names[name] = 1
    return max_key(names)[0]   #max_key returns list of dictionary entries died for the max value
      
def GetLeaves(tree, name):
  nodes = []
  if name == 'REF':
    name = ('EFJ','jgi','ADH')
  elif name == 'KRAUS':
    name = ('KRAUS','kraussiana')   #Include ref seq in monotest
  else:
    name = (name)
  for leaf in tree:
    if any(substring in leaf.name for substring in name):
      nodes.append(leaf)
  return nodes
  
def MonoTest(tree, name):
  nodes = GetLeaves(tree, name)
  ancestor = tree.get_common_ancestor(nodes)
  return len(ancestor) - len(nodes)

def MaxLength(tree, name):
  nodes = GetLeaves(tree, name)
  max_length = 0
  for x in range(len(nodes) - 1):
    for y in range(x+1, len(nodes)):
      dist = tree.get_distance(nodes[x], nodes[y])
      if dist > max_length:
        max_length = dist
      y = y + 1
  return max_length

def GetSisters(node):
  clusters = []
  for c in node.get_children():
    taxa = []
    for l in c.get_leaf_names():
      taxa.append(l)
    clusters.append(taxa)
  return clusters
  

if __name__ == "__main__":
   main(sys.argv[1:])


