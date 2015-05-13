#!/Users/HeathOBrien/anaconda/bin/python


import sys, getopt, requests

def main(argv):
  usage = 'Ensembl.py -q <query> -o <outfile> -f [tree | aa | cds | ortho | ortho_dna] [-v]'
  query = ''
  outfile = ''
  outformat = ''
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
      ext = "/genetree/member/id/%s?nh_format=simple" % query
  elif outformat == 'aa':
      ext = "/sequence/id/%s?content-type=text/x-fasta;type=protein" % query
  elif outformat == 'cds':
      ext = "/sequence/id/%s?content-type=text/x-fasta;type=cds" % query
  elif outformat == 'orthologs' or outformat == 'ortho_dna':
      ext = "/homology/id/%s?content-type=application/json;type=orthologues" % query
  else:
      sys.exit("format %s not recognised\n\n%s" % (outformat, usage))
  if verbose_level > 0:
      sys.stderr.write(server+ext+"\n")
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
              ext = "/sequence/id/%s?multiple_sequences=1;content-type=text/x-fasta;type=cds" % query
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
      out_fh.write(r.text + '\n')

def run_query(request):
    r = requests.get(request)
    if not r.ok:
        r.raise_for_status()
        sys.exit("query %s not recognised" % server+ext)
    return r  

def clean_name(name):
    name = name.replace('>', '')
    if len(name.split('_')) > 2:
        name = '_'.join(name.split('_')[:-2])
    elif len(name.split('_')) == 2:
        sys.exit("Name %s contains one underscore" % name)
    if name[-2] == '.' and name[-1].isdigit(): #this will fail if there are more than 9 transcripts, but I think that's OK
        name = name[:-2]
    return name
    
if __name__ == "__main__":
   main(sys.argv[1:])


