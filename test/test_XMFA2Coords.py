from XMFA2Coords import parse_id, parse_coords
import tempfile, os, warnings
from Bio.Alphabet import generic_dna, IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
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

def test_XMFA2Coords():
    """test if coordinate parsing works correctly"""
    align1 = MultipleSeqAlignment([
             SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="1:1-20 + Seq1"),
             SeqRecord(Seq("ACT-CTAGCCAG", generic_dna), id="2:2-21 + Seq2"),
             SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="1:1-20 - Seq3"),
         ])
    coords = parse_coords(align1)
    print len(coords)
    assert len(coords) == 12
    print len(coords[0])
    assert len(coords[0]) == 3
    print coords[0][0]
    assert coords[0][0] == "Seq1:1"
    print coords[3][1]
    assert coords[3][1] == "-"
    print coords[4][1]
    assert coords[4][1] == "Seq2:5"
    print coords[-1][2]
    assert coords[-1][2] == "Seq3:9"
    
