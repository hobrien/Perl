#!/opt/local/bin/python


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
