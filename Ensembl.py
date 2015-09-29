#!/Users/HeathOBrien/anaconda/bin/python


import sys, getopt, requests, warnings, re

def main(argv):
  usage = 'Ensembl.py -q <query> -o <outfile> -f [tree | aa | cds | ortho | ortho_dna] [-v]'
  query = ''
  outfile = ''
  outformat = ''
  global verbose_level
  verbose_level = 0
  try:
     opts, args = getopt.getopt(argv,"hvq:o:f:",["query=", "outfile=", "outformat=", "verbose="])
  except getopt.GetoptError:
     print usage
     sys.exit(2)
  for opt, arg in opts:
     if opt == '-h':
        print usage
        sys.exit()
     elif opt in ("-q", "--query"):
        query = arg
     elif opt in ("-o", "--outfile"):
        outfile = arg
     elif opt in ("-f", "--outformat"):
        outformat = arg
     elif opt in ("-v", "--verbose"):
        verbose_level = 1
  if outformat == 'ortho':
      outformat = 'orthologs'
  server = "http://rest.ensemblgenomes.org"
  
  query = clean_name(query)
  if outformat == 'tree':
      ext = "/genetree/member/id/%s?content-type=text/x-nh;nh_format=simple" % query.split('.')[0]
  elif outformat == 'aa':
      ext = "/sequence/id/%s?content-type=text/x-fasta;type=protein" % query
  elif outformat == 'cds':
      ext = "/overlap/id/%s?feature=cds;content-type=application/json" % query
      r = run_query(server+ext)
      transcriptID = r.json()[0]['id']
      
      #I'm going to introduce some ugly hacks to deal with the fact that a lot of the
      #transcriptIDs don't show up with this query (if I add object_type=transcript it 
      #throws an error. The only solution that I can think of is to try to guess what the 
      #transcriptID is likely to be:
      transcriptID = transcriptID.replace(r'-PA', r'-TA')
      transcriptID = transcriptID.replace(r'-P', r'')
      transcriptID = re.sub(r'_P(\d\d)$', r'_T\1', transcriptID)
      transcriptID = re.sub(r'(chr\d+)P', r'\1T', transcriptID)
        
      ext = "/sequence/id/%s?content-type=text/x-fasta;type=cds;object_type=transcript" % transcriptID
  elif outformat == 'orthologs' or outformat == 'ortho_dna':
      ext = "/homology/id/%s?content-type=application/json;type=orthologues" % query
  else:
      sys.exit("format %s not recognised\n\n%s" % (outformat, usage))
  r = run_query(server+ext)
  if outfile:
      out_fh = open(outfile, 'w')
  else:
      out_fh = sys.stdout

  if 'ortho' in outformat:
      for ortho in r.json()['data'][0]['homologies']:
          if ('ostreococcus' in ortho['target']['species'] or
              'cyanidioschyzon' in ortho['target']['species'] or
              'chlamydomonas' in ortho['target']['species'] 
              ):
              continue
          if  outformat == 'ortho_dna':
              query = ortho['target']['id']
              ext = "/sequence/id/%s?multiple_sequences=1;content-type=text/x-fasta;type=cds;object_type=transcript" % query
              if verbose_level > 0:
                  sys.stderr.write(server+ext+"\n")
              r2 = run_query(server+ext)
              out_fh.write(">" + r2.text.split(">")[1]) #This prints the first CDS (when there are splice variants)
          else:
              out_fh.write(">%s_%s\n%s\n" % (ortho['target']['id'], 
                                         ortho['target']['species'], 
                                         ortho['target']['align_seq'].replace('-',''))
                      )
  else:
      if outformat == 'cds' and dna_seq(r) != 0:
          sys.exit("Query %s returned seq with %i non-nucleotide characters. Likely AA seq" % (server+ext, dna_seq(r)))
      elif outformat == 'aa' and dna_seq(r) == 0:
          sys.exit("Query %s returned seq with 0 non-nucleotide characters. Likely DNA seq" % server+ext)
      else:
          out_fh.write(r.text + '\n')

def dna_seq(result):
    seq = ''.join(result.text.split('\n')[1:])
    return len(re.findall('[^AaCcGgTtNn]', seq))
    
def run_query(request):
    if verbose_level > 0:
        sys.stderr.write(request+"\n")
    r = requests.get(request)
    if not r.ok:
        warnings.warn("query %s not recognised" % request)
        #r.raise_for_status()
        sys.exit()
    return r  

def clean_name(name):
    name = name.replace('>', '')
    if len(name.split('_')) > 2:
        name = '_'.join(name.split('_')[:-2])
    elif len(name.split('_')) == 2:
        warnings.warn("Name %s contains one underscore" % name)
    return name
    
def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s: %s\n' % (filename, lineno, category.__name__, message)

    
if __name__ == "__main__":
   warnings.formatwarning = warning_on_one_line
   main(sys.argv[1:])


