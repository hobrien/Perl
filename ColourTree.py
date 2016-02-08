#!/Users/HeathOBrien/anaconda/bin/python

"""
This will colour-code taxa on trees according to the Selaginella species
"""

import sys, getopt, string, warnings
import MySQLdb as mdb
from ete2 import Tree, TreeStyle, TextFace, NodeStyle, ImgFace
from os import path  
from Heathpy import get_colour

def main(argv):
  treefilename = ''
  outfilename = ''
  database = 'test'
  usage = 'ColourTree.py -t <treefile> -o <outfile>'
  try:
    opts, args = getopt.getopt(argv,"ht:o:d:",["tree=","out=","db="])
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
    elif opt in ("-d", "--db"):
       database = arg

  tree = Tree(treefilename)
  
  if not outfilename:
    outfilename = treefilename.replace(".nwk", ".pdf")
  con = mdb.connect('localhost', 'root', '', database);
  with con:
    cur = con.cursor()
    for leaf in tree:
      name = leaf.name
      if name.find('KRAUS') > -1:
        color = 'Green'  
      elif name.find('MOEL') > -1:
        color = 'Red'  
      elif name.find('UNC') > -1:
        color = 'Orange'  
      elif name.find('WILD') > -1:
        color = 'MediumBlue'  
      else:
        cur.execute("Select Species.Genus, Species.Species FROM Species, Sequences WHERE Species.abbreviation = Sequences.species AND Sequences.seqid = %s", name)
        try:
          (genus, species) = cur.fetchone()
          leaf.name = '_'.join((genus[0], species, leaf.name))
        except TypeError, e:
          print e
          continue
        if leaf.name.find('kraussiana') > -1:
          color = 'LightGreen'  
        elif leaf.name.find('willdenowii') > -1:
          color = 'SteelBlue' 
        else:
          color = 'Black'
      label = TextFace(leaf.name, fgcolor=color, fsize=16)
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
          expression_label= TextFace(' %s ' % normalized[x], fsize=16)
          expression_label.background.color = get_colour(vsd[x])
          expression_label.border.width = 1
          expression_label.margin_left, expression_label.margin_right, expression_label.margin_top, expression_label.margin_bottom = 1,1,2,1
          # This isn't working right : ( expression_label.border.width=1
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
    #ts.scale = 1500
    ts.show_leaf_name = False
    #ts.show_branch_support = True
    leg_file = path.join(path.expanduser('~'), 'Perl', 'Modules', 'TreeLegend.png')   
    leg_face= ImgFace(img_file=leg_file)
    leg_face.margin_left, leg_face.margin_top = 5, 5
    ts.legend.add_face(leg_face, column=1)
    ts.legend_position=1

    title_face = TextFace(text=file.split('.')[0])
    title_face.margin_left, title_face.margin_top = 10, 5
    ts.title.add_face(title_face, column=1)
    (ts.margin_left, ts.margin_right) = (5,5)
    tree.render(file, tree_style=ts, w=6000, units='mm')
    #tree.show(tree_style=ts)
          
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
    if node.support < 0.9 or node.is_leaf() or node.is_root():
      node.set_style(non_sig)
    else:
      node.set_style(sig)

if __name__ == "__main__":
   main(sys.argv[1:])
