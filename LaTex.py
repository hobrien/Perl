#!/usr/local/bin/python
# -*- coding: utf-8 -*-


import sys, re, subprocess
from os import system
from os.path import expanduser

def ConvertCiteKeys(line):
  """This regex will combine adjacent Papers citekeys, 
  ie {obrien:2010ab}{obrien:2012cd} -> {obrien:2010ab,obrien:2012cd}"""
  line = re.sub(r'([A-Z][a-z]+:\d{4}[a-z]{2})\}\{([A-Z][a-z]+:\d{4}[a-z]{2})', r'\1\2', line)
  
  """This regex will add \citep in front of Papers citekeys (It will fail if suffixes or 
  prefixes are included. I need to add functionality for this)"""
  line = re.sub(r'((\\cite)?\{(([A-Z][a-z]+:\d{4}[a-z]{2},? ?)+)\})', r'\citep\1', line)

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
  			 'Phaeodactylum tricornutum', 'Pseudo-nitzchia multiseries', 'Seminavis robusta', 
  			 'Thalassiosira pseudonana', 'Thalassiosira punctigera',
             'Zonaria tournefortii', 'Z. tournefortii', 'Zonaria']
  for name in species:
    line = re.sub(r'\s%s\s' % name, " \\\\textit{%s} " % name, line)
  return line

def ConvertSymbols(line):
  line = line.replace('µ', "\\textmu ")
  return line

#TODO: add absolute path names for Papers.bib and journal-of-evolutionary-biology.csl
#      option to specify output style

def AddURL(line):
  """ adds \url tags to urls (unless already present). It searches for 'http://' followed 
  a string of valid url characters (including '%xx' unicode triplets) and ending with a
  space or parenthesis/bracket (allowing parens within urls but not at the end
  """
  line = re.sub(r'(\\url{)?(http://([!#$&-;=?-[\]_a-z~]|%[0-9a-fA-F]{2})+)}?[^() [\]]', r'\url{\2}', line)
  return line
  

def pandoc(infile, outfile):
    latex_dir = expanduser('~/Documents/LaTex/')
    print latex_dir
    # TODO manage pandoc errors, for example exit status 43 when citations include Snigowski et al. 2000
    options = ["pandoc", "--bibliography=" + latex_dir + "Papers.bib", "--csl=" + latex_dir + "journal-of-evolutionary-biology.csl", "--output=" + outfile, infile]
    return subprocess.check_call(options)
    #system("pandoc --parse-raw --bibliography=Papers.bib --csl=journal-of-evolutionary-biology.csl --to=latex %s >> %s" % (infile, outfile))
    
#TODO move header to a separate file

header = """\documentclass[a4paper, 12pt, oneside]{article}   	% use "amsart" instead of "article" for AMSLaTeX format
\\usepackage{geometry}                		% See geometry.pdf to learn the layout options. There are lots.
\\geometry{letterpaper}                   		% ... or a4paper or a5paper or ... 
%\geometry{landscape}                		% Activate for for rotated page geometry
%\usepackage[parfill]{parskip}    		% Activate to begin paragraphs with an empty line rather than an indent
\\usepackage{graphicx}				% Use pdf, png, jpg, or epsÂ§ with pdflatex; use eps in DVI mode
\\usepackage{amssymb}
\\usepackage{textgreek}
\usepackage{hyperref}
\\title{Cystoseira}
\\author{Heath E. O'Brien\\\\
        School of Biological Sciences, University of Bristol\\\\
        Woodland Road, Bristol BS8 1UG\\\\
        \\texttt{heath.obrien@gmail.com}}
%\date{}							% Activate to display a given date or no date

\\begin{document}
\\maketitle

"""

#TODO: add options to specify file names for infile and outfile

"""Due to SERIOUS limitations in Pandoc (at least how I'm using it), this needs to be done
in 2 passes. One with Pandoc to format the references, then one with pdflatex to do all 
other formatting.

Not only does pandoc ignore the LaTex control statements, it tends to garble them, so I am
moving as many of them as possible to the second pass"""
in_fh = open(sys.argv[1], 'r')
temp_fh = open('temp.tex', 'w')

temp_fh.write('\\begin{document}\n')
for line in in_fh:
  line = line.strip()
  line = ConvertCiteKeys(line)
  #print line
  temp_fh.write(line + '\n')    
in_fh.close()
temp_fh.write("\n\section*{References}\n")   
temp_fh.write("\end{document}\n")
temp_fh.close()           
pandoc("temp.tex", "pandoc.tex")

pandoc_fh = open("pandoc.tex", 'r')
out_fh = open("out.tex", 'w')
out_fh.write(header)
for line in pandoc_fh:
  line = line.strip()
  line = AddItalics(line)
  line = ConvertSymbols(line)
  line = AddURL(line)
  out_fh.write(line + '\n')
pandoc_fh.close()
out_fh.write("\n\end{document}\n")
out_fh.close()
#system("rm temp.tex pandoc.tex")
system("pdflatex out.tex")
