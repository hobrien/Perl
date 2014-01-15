#!/usr/local/bin/python


import sys, getopt
from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Nexus import Nexus
from Heathpy import write_phylip

def main(argv):
  usage = 'ConvertAln -i <infile> -x <informat> -o <outfile> -f <outformat>'
  infile = ''
  informat = ''
  outfile = ''
  outformat = ''
  try:
     opts, args = getopt.getopt(argv,"hi:x:o:f:",["infile=", "informat=", "outfile=", "outformat="])
  except getopt.GetoptError:
     print usage
     sys.exit(2)
  for opt, arg in opts:
     if opt == '-h':
        print usage
        sys.exit()
     elif opt in ("-i", "--infile"):
        infile = arg
     elif opt in ("-x", "--informat"):
        informat = arg
     elif opt in ("-o", "--outfile"):
        outfile = arg
     elif opt in ("-f", "--outformat"):
        outformat = arg
  if not infile:
    sys.exit("must specify infile! %s" % usage)

  if not outformat:
    sys.exit("must specify format to convert to! %s" % usage)
  
  if not informat:
    informat = guess_format(infile)

  if informat == 'phylip':
    informat = 'phypli-relaxed'  
  print outfile
  if not outfile:
    if '.' in infile:
      outfile = '.'.join((infile.split('.')[:-1] + [get_extension(outformat)]))
    else:
      outfile = '.'.join((infile, get_extension(outformat)))
          
  if outformat == 'nexus':
    alignment=AlignIO.read(infile, informat, alphabet=IUPAC.ambiguous_dna)
    write_nexus(alignment, outfile)
  
  elif outformat == 'phylip':
    alignment=AlignIO.read(infile, informat, alphabet=IUPAC.ambiguous_dna)
    out_fh = open(outfile, 'w')
    write_phylip(alignment, out_fh)
    out_fh.close()
  
  else:
    AlignIO.convert(infile, informat, outfile, outformat, alphabet=IUPAC.ambiguous_dna)

def write_nexus(alignment, outfile):
  minimal_record = "#NEXUS\nbegin data; dimensions ntax=0 nchar=0; format datatype=%s; end;" % "dna"
  n = Nexus.Nexus(minimal_record)
  n.alphabet = alignment._alphabet
  for record in alignment:
    n.add_sequence(record.id, record.seq.tostring())
    n.write_nexus_data(outfile, interleave=False)

def guess_format(file):
  ext = file.split('.')[-1]
  if ext[0] == 'f' and ext[-1] == 'q': format = "fastq"
  elif ext[0:2] == 'fa': format = "fasta"
  elif ext == 'xmfa': format = "xmfa"
  elif ext == 'phy': format = "phylip-relaxed"
  elif ext == 'aln': format = "clustal"
  elif ext == 'gb': format =  "genbank"
  elif ext == 'nxs': format =  "nexus"
  elif ext == 'nex': format =  "nexus"
  elif ext == 'maf': format =  "maf"
  else: sys.exit("extension not recognized! Use -x option to specify explicitly")
  return format

def get_extension(format):
  if format == "fasta": extension = "fa"
  elif format == "xmfa": extension = "xmfa"
  elif format == "fastq": extension = "fq"
  elif format == "phylip": extension = "phy"
  elif format == "clustal": extension = "aln"
  elif format == "genbank": extension = "gbk"
  elif format == "nexus": extension = "nex"  
  else: sys.exit("format not recognized! Please select one of 'clustalw', 'fasta', 'fastq', 'genbank', 'nexus', 'phylip' or 'xmfa'")
  return extension
  
if __name__ == "__main__":
  if sys.argv[1] == 'self-test':
    from os import path, system
    print "Converting fasta to phylip..."
    infile = path.join(path.expanduser("~"), 'Perl', 't', 'ConvertAln.t.fa')
    main(('-i', infile, '-f', 'phylip'))
    system('open %s' % path.join(path.expanduser("~"), 'Perl', 't', 'ConvertAln.t.phy'))
    print "Converting phylip to nexus..."
    infile = path.join(path.expanduser("~"), 'Perl', 't', 'ConvertAln.t.phy')
    main(('-i', infile, '-f', 'nexus'))
    system('open %s' % path.join(path.expanduser("~"), 'Perl', 't', 'ConvertAln.t.nex'))
    print "Converting nexus to clustal..."
    infile = path.join(path.expanduser("~"), 'Perl', 't', 'ConvertAln.t.nex')
    main(('-i', infile, '-f', 'clustal'))
    system('open %s' % path.join(path.expanduser("~"), 'Perl', 't', 'ConvertAln.t.aln'))

  
  else:
    main(sys.argv[1:])

