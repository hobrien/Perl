from Heathpy import make_hash
from StringIO import StringIO
import tempfile, os, warnings

test_dir = os.path.join(os.path.expanduser('~'), 'Perl', 'test')
def test_make_hash1():
  data = """seq1 500
  seq2 1500
  seq3 2500"""
  hash = make_hash(StringIO(data))
  assert len(hash) == 3
  assert 'seq1' in hash
  assert 'seq2' in hash
  assert 'seq3' in hash
  assert 'seq4' not in hash
    
def test_make_hash2():
  data = """SeqID length
  seq1 500
  seq2 1500
  seq3 2500"""
  hash = make_hash(StringIO(data),2,1)
  assert len(hash) == 3
  assert '500' in hash
  assert '1500' in hash
  assert '2500' in hash
  assert '2000' not in hash

def test_make_hash3():
  table = os.path.join(test_dir, 'test_table.txt')
  print table
  fh = open(table, 'r')
  hash = make_hash(fh,1,1, '\t')
  assert len(hash) == 3
  assert '691545' in hash
  assert '691883' in hash
  assert '699770' in hash
  assert '2000' not in hash
  fh.close()
  
def test_make_hash4():
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    table = os.path.join(test_dir, 'test_table.csv')
    fh = open(table, 'r')
    hash = make_hash(fh, sep=',')
    fh.close()
    assert len(hash) == 3
    assert '691545' in hash
    assert '691883' in hash
    assert '699770' in hash
    assert '2000' not in hash
    assert len(w) == 1
    assert "multiple entries" in str(w[-1].message)