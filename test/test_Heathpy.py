from Heathpy import make_hash
from StringIO import StringIO
import tempfile, os, warnings

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
