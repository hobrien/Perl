#!/usr/local/bin/python
# -*- coding: utf-8 -*-

"""I am rapidly hitting the limits of what I can do with pandoc/bibtex and I hate having
to maintain two databases. I'm going to remove the citation formatting from this and just
do that with the Papers citation support. This means all I have to do is one pass with
LaTex"""

import sys, re, subprocess
from os import system, path

def ConvertCiteKeys(line):
  """This regex will combine adjacent Papers citekeys, 
  ie {obrien:2010ab}{obrien:2012cd} -> {obrien:2010ab,obrien:2012cd}"""
  line = re.sub(r'([A-Z][a-z]+:\d{4}[a-z]{2})\}\{([A-Z][a-z]+:\d{4}[a-z]{2})', r'\1\2', line)
  
  """This regex will add \citep in front of Papers citekeys (It will fail if suffixes or 
  prefixes are included. I need to add functionality for this)"""
  line = re.sub(r'(\\cite)?(\{(([A-Z][a-z]+:\d{4}[a-z]{2},? ?)+)\})', r'\citep\2', line)

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
    line = re.sub(r'\b%s\b' % name, " \\\\textit{%s} " % name, line)
  return line

def ConvertSymbols(line):
  line = line.replace('µ', "\\textmu ")
  line = line.replace('&', "\& ")
  return line

def AddURL(line):
  """ adds \url tags to urls (unless already present). It searches for 'http://' followed 
  a string of valid url characters (including '%xx' unicode triplets) and ending with a
  space or parenthesis/bracket (allowing parens within urls but not at the end
  """
  line = re.sub(r'(\\url{)?(http://([!#$&-;=?-[\]_a-z~]|%[0-9a-fA-F]{2})+)}?[^() [\]]', r'\url{\2}', line)
  return line
  
  # TODO: deal with urls that don't contain 'http://'
  

def pandoc(infile, outfile):
    latex_dir = path.expanduser('~/Documents/LaTex/')
    print latex_dir
    options = ["pandoc", "--bibliography=" + latex_dir + "Papers.bibtex", "--csl=" + latex_dir + "journal-of-evolutionary-biology.csl", "--output=" + outfile, infile]
    return subprocess.check_call(options)

    # TODO: manage pandoc errors, for example exit status 43 when citations include Snigowski et al. 2000
    # TODO: add option to specify output style
    

"""Due to SERIOUS limitations in Pandoc (at least how I'm using it), this needs to be done
in 2 passes. One with Pandoc to format the references, then one with pdflatex to do all 
other formatting.

Not only does pandoc ignore the LaTex control statements, it tends to garble them, so I am
moving as many of them as possible to the second pass"""
infile = sys.argv[1]
basename, ext = path.splitext(infile)
if ext == 'tex':
  outfile = basename + '_out.tex'
else:
  outfile = basename + '.tex'
  
in_fh = open(infile, 'r')

#read first line of infile to obtain title
title = in_fh.readline()

header = """\documentclass[a4paper, 12pt, oneside]{article}   	% use "amsart" instead of "article" for AMSLaTeX format
\\usepackage[utf8]{inputenc}
\\usepackage{geometry}                		% See geometry.pdf to learn the layout options. There are lots.
\\geometry{letterpaper}                   		% ... or a4paper or a5paper or ... 
%\geometry{landscape}                		% Activate for for rotated page geometry
%\usepackage[parfill]{parskip}    		% Activate to begin paragraphs with an empty line rather than an indent
\\usepackage{graphicx}				% Use pdf, png, jpg, or epsÂ§ with pdflatex; use eps in DVI mode
\\usepackage{amssymb}
\\usepackage{textgreek}
\usepackage{hyperref}
\\title{""" + title + """}
\\author{Heath E. O'Brien\\\\
        School of Biological Sciences, University of Bristol\\\\
        Woodland Road, Bristol BS8 1UG\\\\
        \\texttt{heath.obrien@gmail.com}}
%\date{}							% Activate to display a given date or no date

\\begin{document}
\\maketitle

"""
#TODO move header to a separate file


out_fh = open(outfile, 'w')
out_fh.write(header)
for line in in_fh:
  line = line.strip()
  line = AddItalics(line)
  line = ConvertSymbols(line)
  line = AddURL(line)
  out_fh.write(line + '\n')
in_fh.close()
out_fh.write("\n\end{document}\n")
out_fh.close()
system("pdflatex %s" % outfile)
#system("rm temp.tex pandoc.tex %s.tex %s.log %s.aux %s.out" % (basename, basename, basename, basename))
