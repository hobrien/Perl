#!/opt/local/bin/python
import sys
from Heathpy import parse_xmfa

def parse_id(id):
  id_parts = id.split(' ')
  if len(id_parts) != 3:
    raise TypeError("id format should be '#:left-right +/- name")
  name = id_parts[2]
  strand = id_parts[1]
  coords = id_parts[0].split(':')[1].split('-')
  if len(coords) != 2:
    raise TypeError("id format should be '#:left-right +/- name")
  if strand == '+':
    start = int(coords[0])
  elif strand == '-':
    start = int(coords[1])
  else:
    raise TypeError("strand information should be either '+' or '-'")
  return (start, strand, name)

def parse_coords(aln):
  """This takes an alignment where sequence ids are formatted using Mauve formatting
  and produces a list of matching coordinates for each sequence, with the name of the
  sequence separated from the position with a colon."""
  coord_array = []    # list of homologous positions
  column_coords = []  # current coordinate position for each seq in alignment
  strands = [] # list of strand information for each seq in alignment
  names = [] # list of sequence names
  for seq in aln:
    (start, strand, name) = parse_id(seq.id)
    names.append(name)
    strands.append(strand)
    column_coords.append(start)
  for x in range(len(aln[0])):
    row_coords = []
    slice = aln[:, x]
    for y in range(len(slice)):
      if slice[y] != '-':
        row_coords.append(names[y] + ':' + str(column_coords[y]))
        if strands[y] == '+':
          column_coords[y] += 1
        else:
          column_coords[y] -= 1
      else:
        row_coords.append('-')  
    coord_array.append(row_coords)
  return coord_array

if __name__ == "__main__":
   fh = open(sys.argv[1], 'r')
   for alignment in parse_xmfa(fh):
     if len(alignment) < 2:
       continue
     for row in parse_coords(alignment):
       print ", ".join(row)
