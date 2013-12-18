#!/usr/local/bin/python
# -*- coding: utf-8 -*-


import sys, re, subprocess

def ConvertCiteKeys(line):
  """This regex will combine adjacent Papers citekeys, 
  ie {obrien:2010ab}{obrien:2012cd} -> {obrien:2010ab,obrien:2012cd}"""
  line = re.sub(r'([A-Z][a-z]+:\d{4}[a-z]{2})\}\{([A-Z][a-z]+:\d{4}[a-z]{2})', r'\1\2', line)
  
  """This regex will add \citep in front of Papers citekeys (It will fail if suffixes or 
  prefixes are included. I need to add functionality for this)"""
  line = re.sub(r'(\{(([A-Z][a-z]+:\d{4}[a-z]{2},? ?)+)\})', r'\citep\1', line)

  """This regex is to deal with author suppression (ie OBrien et al. (2012) showed that...
  It will fail if the citekey has multiple papers, though I can't imagine a situation where
  this would happen
  
  \citet would be a better way to handle this situation because it will handle things like 
  formatting et al for the correct style automaticaly, but it would be harder to implement"""
  line = re.sub(r'(\{\*)([A-Z][a-z]+:\d{4}[a-z]{2}\})', r'\citeyear{\2', line)
  
  """This regex will handle year suppression. See above for limitations"""
  line = re.sub(r'(\{[A-Z][a-z]+:\d{4}[a-z]{2})(\*\})', r'\citeauthor\1}', line)
  return line


"""This works OK, though there may be a more efficient way to do it using regex"""
def AddItalics(line):
  species = ['Cysoseira stricta', 'C. stricta', 'Cysoseira', 'Dictyota dichotoma', 'D. dichotoma', 'Dictyota', 
           'Zonaria tournefortii', 'Z. tournefortii', 'Zonaria']
  for name in species:
    line = re.sub(r'\s%s\s' % name, " \\\\textit{%s} " % name, line)
  return line


"""This is working reasonably well, except for the fact that almost all of the formatting
info in the header is ignored by pandoc. I haven't been able to figure out why this is,
but it will need to be solved for any of this to be useful

TODO: add absolute path names for Papers.bib and journal-of-evolutionary-biology.csl
      option to specify output style
"""

def pandoc(infile, outfile):
    # TODO manage pandoc errors, for example exit status 43 when citations include Snigowski et al. 2000
    options = ["pandoc", "--bibliography=Papers.bib", "--csl=journal-of-evolutionary-biology.csl", "-o", outfile, infile]
    return subprocess.check_call(options)
    

header = """\documentclass[11pt, oneside]{article}   	% use "amsart" instead of "article" for AMSLaTeX format
\\usepackage{geometry}                		% See geometry.pdf to learn the layout options. There are lots.
\\geometry{letterpaper}                   		% ... or a4paper or a5paper or ... 
%\geometry{landscape}                		% Activate for for rotated page geometry
%\usepackage[parfill]{parskip}    		% Activate to begin paragraphs with an empty line rather than an indent
\\usepackage{graphicx}				% Use pdf, png, jpg, or epsÂ§ with pdflatex; use eps in DVI mode
\\usepackage{amssymb}
\\title{Cystoseira}
\\author{Heath E. O'Brien}
%\date{}							% Activate to display a given date or no date

\\begin{document}
\\maketitle

%\section{}
%\subsection{}
"""

"""TODO: add options to specify file names for infil and outfile"""
in_fh = open(sys.argv[1], 'r')
out_fh = open('out.tex', 'w')

out_fh.write(header)
for line in in_fh:
  line = line.strip()
  line = ConvertCiteKeys(line)
  line = AddItalics(line)
  #print line
  out_fh.write(line + '\n')     
out_fh.write("\end{document}")
out_fh.close()           

pandoc("out.tex", "out.pdf")