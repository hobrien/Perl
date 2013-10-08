#!/usr/local/bin/python

import sys, getopt
from ete2 import Tree
import re
from os import listdir

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'CheckMonoplyly.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'CheckMonoplyly.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputdir = arg
#      elif opt in ("-o", "--ofile"):
#         outputfile = arg
   species_list = ['UNC', 'MOEL', 'WILD', 'KRAUS', 'REF']
   print "Tree\tGenus_mono\tGenus_max_legnth\tNum_mono\tSpecies_max_length\tCongruence",
   for species in species_list:
     print "\t%s_mono\t%s_length" % (species, species),
   print ""  
   
   for inputfile in listdir(inputdir):
     if re.search('nwk', inputfile):
       input = inputdir + "/" + inputfile
       t = Tree(input)
       R = t.get_midpoint_outgroup()
       t.set_outgroup(R)
       output = []
       global_max = 0;
       num_mono = 0;
       for species in species_list:
         leaves = GetLeaves(t, species)
         if len(leaves) < 2:
           num_mono += 1
           output += [ 'mono', 0]
         else:
           mono = MonoTest(t, species)
           if mono == 0:
             mono = "mono"
             num_mono += 1
           length = MaxLength(t, species)       
           output += [ mono, length]
           if length > global_max:
             global_max = length
       congruence = 'incongruent'
       if MonoTest(t, '(UNC)|(WILD)') == 0 and MonoTest(t, '(UNC)|(WILD)|(MOEL)') == 0:
         congruence = 'congruent'
       genus_mono = 'poly'
       if MonoTest(t, '(UNC)|(WILD)|(MOEL)|(EFJ)|(jgi)|(ADH)|(KRAUS)') == 0:
         genus_mono = 'mono'
       genus_length = MaxLength(t, '(UNC)|(WILD)|(MOEL)|(EFJ)|(jgi)|(ADH)|(KRAUS)')
       
       output = [ inputfile, genus_mono, genus_length, num_mono, global_max, congruence] + output
       print '\t'.join(str(x) for x in output)


def GetLeaves(tree, name):
  nodes = []
  if name == 'REF':
     name = '(EFJ)|(jgi)|(ADH)'
  for leaf in tree:
    if re.search(name, leaf.name):
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


if __name__ == "__main__":
   main(sys.argv[1:])


