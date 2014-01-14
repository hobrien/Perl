#!/usr/local/bin/python

import subprocess

def get_taxid(id):
  """takes GI number, Genbank accession number or full genbank id and looks up taxID number
  using blastdbcmd"""
  options = ("blastdbcmd -entry '" + id + "' -db nt -outfmt %T")
  proc = subprocess.Popen(options, stdout=subprocess.PIPE, shell=True)
  staxids = proc.communicate()
  return int(staxids[0].strip())
  


