#!/opt/local/bin/python
"""

MakeConsesus.py 5 Sept 2013

SYNOPSIS

MakeConsensus.py -t <treefile> -a <alignmentfile> -s <species> -c <consensusfile> -l <logfile>

Options:
 -t treefile
 -a alignment file (relaxed Phylip format)
 -s species (KRAUS, UNC, WILD, MOEL)
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
(in case of a tie, one is picked arbitrarily). If multiple consensus sequences are derived
from the same comp, numbers are appended to the names (ie comp1_cons, comp1_1_cons, 
com1_2_cons, etc). This works across multiple runs of the script, as long as they are
launched from the same folder

AUTHOR

 Heath E. O'Brien (heath.obrien-at-gmail-dot-com)
"""

import sys  
import getopt 
from os import path, system
from ete2 import Tree  
import cPickle as pickle
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
   species_name = ''
   try:
      opts, args = getopt.getopt(argv,"ht:a:c:l:s:",["tree=","alignment=","species=", "consensus=","log="])
   except getopt.GetoptError:
      print 'MakeConsensus.py -t <treefile> -a <alignmentfile> -s <species> -c <consensusfile> -l <logfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'MakeConsensus.py -t <treefile> -a <alignmentfile> -s <species> -c <consensusfile> -l <logfile>'
         sys.exit()
      elif opt in ("-t", "--tree"):
         treefile = arg
      elif opt in ("-a", "--alignment"):
         alnfile = arg
      elif opt in ("-c", "--consensus"):
         consfilename = arg
      elif opt in ("-l", "--log"):
         logfilename = arg
      elif opt in ("-s", "--species"):
         species_name = arg
   
   
   #make a list of names to ensure that there are no duplicates.
   #if saved list exists, use it instead
   namefile = "names.p"
   if path.exists(namefile):
     name_list = pickle.load( open( namefile, "rb" ) )
   else:
     name_list = {}
   name_num = 1
   
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
   
   #Add "species" feature to leaves on tree for easy lookup
   t = add_species(t)
   
   #cycle through monophyletic groups and build consensus sequences
   for node in t.get_monophyletic(values=[species_name], target_attr="species"):
     species_seqs = []
     for leaf in node:
       if 'kraussiana' not in leaf.name:                   #exclude reference sequence from consensus
         species_seqs.append(aln[aln_index[leaf.name]])
     if len(species_seqs) > 1:                             #create consensus sequence from sequence list
       species_aln = MultipleSeqAlignment(species_seqs)
       consensus = smart_consensus(MultipleSeqAlignment(pad_ends(species_aln))) #use pad_ends to convert end gaps to ambiguous
       consensus = pad_ends([consensus], '-', 'N')[0]      #convert back to gap characters (need to convert seq to a list)
       name = pick_name(species_seqs) + "_cons"
       while name in name_list:
         name = name.split("_")[0] + "_" + str(name_num) + "_cons"
         name_num = name_num + 1
       name_list[name] = 1
       name_num = 1
       consensus =  SeqRecord(Seq(str(consensus).replace('-','')),
                              id=name, 
                              description=name)
     else:                                                #for singletons, just write ungapped version of original
       consensus = species_seqs[0]
       consensus.seq = consensus.seq.ungap("-")                    
       name = consensus.id
       
     #write consensus sequence to consensus file
     consfile = open(consfilename, "a")
     consfile.write(consensus.format("fasta"))
     consfile.close()
       
     #write summary of seqs represented by cons to log file
     cluster_num = path.split(alnfile)[1].split(".")[0] #get cluster number
     logfile = open(logfilename, "a")
     for seq in species_seqs:
       logfile.write("%s\n" % ", ".join([seq.id, name, cluster_num]))
     logfile.close()

   #save list of names for future uses of script
   pickle.dump( name_list, open( namefile, "wb" ) )

def add_species(tree): 
   species_list = ['UNC', 'MOEL', 'WILD', 'KRAUS']
   for leaf in tree:
     to_add = "other"
     for species in species_list:
       if species in leaf.name:
         to_add = species
     if "kraussiana" in leaf.name:       #Include reference sequences from 1kp project
       to_add = "KRAUS"
     leaf.add_features(species=to_add)
   return tree
   
def pick_name(seqs):
  names = {}
  for seq in seqs:
    name = seq.id.split("_")[0]
    if name in names:
      names[name] = names[name] + 1
    else:
      names[name] = 1
  return max_key(names)[0]   #max_key returns list of dictionary entries tied for the max value
      

if __name__ == "__main__":
   main(sys.argv[1:])


