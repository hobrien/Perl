from XMFA2Coords import parse_id
import tempfile, os, warnings
from nose.tools import *
from nose import with_setup

"""This is to be used with nosetests: 'nosetests ~/Perl' """
test_dir = os.path.join(os.path.expanduser('~'), 'Perl', 'test')

def test_parse_id():
  (start, strand, name) = parse_id('1:1-20 + /Users/HeathOBrien/Bioinformatics/Psuedomonas/AE016853.gbk')
  assert start == 1
  assert strand == '+'
  assert name == '/Users/HeathOBrien/Bioinformatics/Psuedomonas/AE016853.gbk'
  (start, strand, name) = parse_id('1:1-20 - /Users/HeathOBrien/Bioinformatics/Psuedomonas/AE016853.gbk')
  assert start == 20
  assert strand == '-'

@raises(TypeError)  
def test_parse_id2():
  """test if error handling works correctly when id not formatted correctly"""
  (start, strand, name) = parse_id('1:1-20,1,/Users/HeathOBrien/Bioinformatics/Psuedomonas/AE016853.gbk')

@raises(TypeError)  
def test_parse_id3():
  """test if error handling works correctly when strand info not formatted correctly"""
  (start, strand, name) = parse_id('1:1-20 1 /Users/HeathOBrien/Bioinformatics/Psuedomonas/AE016853.gbk')
