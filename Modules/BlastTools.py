#!/usr/local/bin/python

import subprocess

def get_taxid(id):
  #id = "'%s'" % id
  options = ("blastdbcmd -entry '" + id + "' -db nt -outfmt %T")
  proc = subprocess.Popen(options, stdout=subprocess.PIPE, shell=True)
  staxids = proc.communicate()
  return staxids[0].strip()
  


