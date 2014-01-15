from Heathpy import make_hash, remove_dots
from StringIO import StringIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import tempfile, os, warnings
from nose.tools import *

test_dir = os.path.join(os.path.expanduser('~'), 'Perl', 'test')

def test_make_hash1():
  """test if able to skip header line"""
  table = os.path.join(test_dir, 'test_table.txt')
  print table
  fh = open(table, 'r')
  hash = make_hash(table,1,1, '\t')
  assert len(hash) == 3
  assert '691545' in hash
  assert '691883' in hash
  assert '699770' in hash
  assert '2000' not in hash
  fh.close()
  
def test_make_hash2():
  """test for warning when there are multiple entires for an id"""
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    table = os.path.join(test_dir, 'test_table.csv')
    fh = open(table, 'r')
    hash = make_hash(table, sep=',')
    fh.close()
    assert len(w) == 1
    assert "multiple entries" in str(w[-1].message)
    
def test_make_hash3():
  """test if sniffer able to identify delimiter"""
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    table = os.path.join(test_dir, 'test_table.csv')
    fh = open(table, 'r')
    hash = make_hash(table)
    fh.close()
    assert len(hash) == 3
    assert '691545' in hash
    assert '691883' in hash
    assert '699770' in hash
    assert '2000' not in hash

def test_remove_dots():
  align1 = MultipleSeqAlignment([
             SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="Alpha"),
             SeqRecord(Seq("ACT-CTAGC.AG", generic_dna), id="Beta"),
             SeqRecord(Seq("ACTGCTAGDTAG", generic_dna), id="Gamma"),
         ])
  assert str(remove_dots(align1)[1].seq) == "ACT-CTAGCTAG"

@raises(AssertionError)  
def test_remove_dots2():
  """test if error handling works correctly when dot in reference seq"""
  align1 = MultipleSeqAlignment([
             SeqRecord(Seq("A.TGCTAGCTAG", generic_dna), id="Alpha"),
             SeqRecord(Seq("ACT-CTAGC.AG", generic_dna), id="Beta"),
             SeqRecord(Seq("ACTGCTAGDTAG", generic_dna), id="Gamma"),
         ])
  remove_dots(align1)