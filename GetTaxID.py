#!/usr/local/bin/python

import fileinput
from BlastTools import get_taxid

  
for line in fileinput.input():
  if "Hit_ID" in line: continue
  fields = line.split()
  taxid =  get_taxid(fields[2])
  if taxid > 0:
    print '\t'.join(fields[0], taxid)

