#!/usr/local/bin/python
"""

GetClades.py 12 Sept 2013

SYNOPSIS

MakeClades.py -t <treefile> -o <outfile> <species list>

Options:
 -t treefile
 -o outfile for list of clades for each species
 <species list>

DESCRIPTION
Takes a tree file and a list of species, finds all clades exclusive to these species
and produces a list of which clades each species occurs in

NOTES
results are APPENDED to the outfile so that data can be collected form multiple trees

AUTHOR

 Heath E. O'Brien (heath.obrien-at-gmail-dot-com)
"""

import sys  
import getopt 
from os import path, system
from ete2 import Tree  

def main(argv):
   treefile = ''
   species_list = argv[2:]
   try:
      opts, args = getopt.getopt(argv[0:2],"ht:o:",["tree=","alignment=","species=", "consensus=","log="])
   except getopt.GetoptError:
      print 'GetClades.py -t <treefile> -o <outfile> <species>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'GetClades.py -t <treefile> -o <outfile> <species>'
         sys.exit()
      elif opt in ("-t", "--tree"):
         treefile = arg
   
   #create outfilename from treefile
   outfilename = treefile.replace(".nwk", "_clades.txt")
   

   #skip tree files with a size of zero (usually caused by too few sequences to make a tree)
   if path.getsize(treefile) == 0:
     sys.exit()      


     
   #read in tree, root by midpoint, and add info about target species
   t = Tree(treefile)
   if len(t) > 3:                 
     R = t.get_midpoint_outgroup()  #midpoint root tree unless there are only three taxa, which was causing an error
     t.set_outgroup(R)
   t = add_species([t] + species_list)
   
   #extract cluster name from treefile (remove path and extension)
   cluster_name = path.split(treefile)[1].split(".")[0] + "_"
   #make dictionary with clade numbers as keys and arrays of species as values
   clade_dict = {}
   clade_num = 0
   for node in t.get_monophyletic(values=["ingroup"], target_attr="species"):
     clade_num = clade_num +1
     clade_dict[clade_num] = [""] * len(species_list)  #create array of empty values with same length as species_list
     for leaf in node:
       for i in range(len(species_list)):
         if species_list[i] in leaf.name:
           clade_dict[clade_num][i] = cluster_name + str(clade_num)

   #open outfile and write results
   outfile_handle = open(outfilename, "w")

   for clade in clade_dict.values():
     outfile_handle.write("\t".join(clade) + "\n")
   outfile_handle.close()
  
      
def add_species(species_list):
  tree = species_list.pop(0)
  for leaf in tree:
    to_add = 'outgroup'
    for species in species_list:
      if species in leaf.name:
        to_add = 'ingroup'
    leaf.add_features(species=to_add)
  return tree

if __name__ == "__main__":
   main(sys.argv[1:])


