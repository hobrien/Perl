#!/usr/local/bin/python

"""
This will colour-code taxa on trees according to the Selaginella species
"""

import sys, getopt, string, warnings
import MySQLdb as mdb
from ete2 import Tree, TreeStyle, TextFace, NodeStyle  
from Heathpy import get_colour

def main(argv):
  treefilename = ''
  outfilename = ''
  usage = 'ColourTree.py -t <treefile> -o <outfile>'
  try:
    opts, args = getopt.getopt(argv,"ht:o:",["tree=","out=",])
    if not opts:
      raise getopt.GetoptError('no opts')
  except getopt.GetoptError:
    print usage
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print usage
       sys.exit()
    elif opt in ("-t", "--tree"):
       treefilename = arg
    elif opt in ("-o", "--out"):
       outfilename = arg

  tree = Tree(treefilename)
  
  if not outfilename:
    outfilename = treefilename.replace(".nwk", ".pdf")
  con = mdb.connect('localhost', 'root', '', 'Selaginella');
  with con:
    cur = con.cursor()
    for leaf in tree:
      name = leaf.name
      if name.find('KRAUS') > -1:
        color = 'Green'  
      elif name.find('kraussiana') > -1:
        color = 'SpringGreen'  
      elif name.find('MOEL') > -1:
        color = 'Red'  
      elif name.find('UNC') > -1:
        color = 'Orange'  
      elif name.find('WILD') > -1:
        color = 'MediumBlue'  
      elif name.find('willdenowii') > -1:
        color = 'SteelBlue'  
      else:
        color = 'Black'
      label = TextFace(leaf.name, fgcolor=color)
      #label.background.color = color
      leaf.name = label
      leaf.add_face(label, column = 0, position="branch-right")
      leaf.add_face(TextFace(' '), column = 1, position="branch-right")
      name = '_'.join(name.split('_')[:-2])
      try:
        cur.execute("SELECT vsd.leaf1, vsd.leaf2, vsd.leaf3, vsd.leaf4 FROM vsd, orfs WHERE vsd.gene_id = orfs.gene_id AND orfs.seqid = %s", name)
        vsd = cur.fetchone()
        cur.execute("SELECT normalized.leaf1, normalized.leaf2, normalized.leaf3, normalized.leaf4 FROM normalized, orfs WHERE normalized.gene_id = orfs.gene_id AND orfs.seqid = %s", name)        
        normalized = cur.fetchone()
        for x in range(4):
          if vsd[x] == 'none':
            continue
          expression_label= TextFace(' %s ' % normalized[x])
          expression_label.background.color = get_colour(vsd[x])   
          leaf.add_face(expression_label, column = x+2, position="branch-right")
      except TypeError:
       continue
    draw_tree(tree, outfilename)   

def draw_tree(tree, file):
    root = tree.get_midpoint_outgroup()
    try:
      tree.set_outgroup(root)
    except:
      pass
    root = tree.get_tree_root()
    root.dist = 0
    add_sig(tree)
    ts = TreeStyle()
    ts.branch_vertical_margin = 1
    ts.scale = 1500
    ts.show_leaf_name = False
    tree.render(file, tree_style=ts, w=3000, units='mm')
    #tree.show(tree_style=ts)
    
def add_faces(leaf, label_info):
      colours = get_colours(label_info)
      y = 0
      for x in range(len(label_info)):
        if x < len(label_info) - 1:
          label_info[x] += ','
        label = TextFace(label_info[x])
        label.margin_left = 5
        label.fgcolor = colours[x]
        if x > 1 and x % 3 == 0:
          y += 3
        leaf.add_face(label, column=x-y+1, position="branch-right")
      
def get_colours(label_info):
  colours = []
  for label in label_info:
    genus = label.split(' ')[0]
    if genus.find('.') != -1:
      colours.append(colours[-1])
    else:
      con = mdb.connect('localhost', 'root', '', 'PhotobiontDiversity');
      with con:
        cur = con.cursor()
        try:
          cur.execute("SELECT phylum FROM Taxonomy WHERE genus= %s", (genus))
          taxon = cur.fetchone()
        except TypeError:
          warnings.warn("No phylum entry for %s" % genus)
        if taxon and taxon[0] == 'Ascomycota':
          try:
            cur.execute("SELECT family FROM Taxonomy WHERE genus= %s", (genus))
            taxon = cur.fetchone()
          except TypeError:
            warnings.warn("No family entry for %s" % genus)
        try:
          cur.execute("SELECT Colour FROM Colours WHERE Taxon= %s", (taxon[0]))
          colour = cur.fetchone()      
          colours.append(colour[0])
        except TypeError:
          warnings.warn("No colour available for %s (%s)" % (genus, taxon))
          colours.append('LightGray')
          
  return colours
    
def add_sig(tree):
  non_sig = NodeStyle()
  non_sig["size"] = 0
  sig = NodeStyle()
  sig["size"] = 5
  sig["fgcolor"] = "black"
  for node in tree.traverse():
    if node.support < 0.9 or node.is_leaf():
      node.set_style(non_sig)
    else:
      node.set_style(sig)

def combine_info(entries):
  host_counts = {}                   #Can include species names of free-living strains
  for (host, species) in entries:
    if host and host != "free-living":
      info = host
    else:
      info = species
    if info in host_counts.keys():
      host_counts[info] += 1
    else:
      host_counts[info] = 1
  out_list = []
  included_genera = []
  names = host_counts.keys()
  names.sort()
  for name in names:
      count = host_counts[name]
      genus = name.split(' ')[0]
      if genus in included_genera:
        name = name.replace(genus,  genus[0] + '.')
      else:
        included_genera.append(genus)
      if count > 1:  
        out_list.append("%s (%s)" % ( name, count))
      else:
        out_list.append(name)

  return out_list
  
if __name__ == "__main__":
   main(sys.argv[1:])
