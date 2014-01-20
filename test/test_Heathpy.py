from Heathpy import make_hash, remove_dots, write_phylip, six_frame_translation, find_CaaX
from StringIO import StringIO
from Bio.Alphabet import generic_dna, IUPAC
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


align1 = MultipleSeqAlignment([
             SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="Alpha"),
             SeqRecord(Seq("ACT-CTAGC.AG", generic_dna), id="Beta"),
             SeqRecord(Seq(".....A......", generic_dna), id="Gamma"),
         ])

def test_remove_dots():
  aln2 = remove_dots(align1)
  print "testing if %s == ACTGCTAGCTAG" % aln2[0].seq       
  assert str(remove_dots(align1)[0].seq) == "ACTGCTAGCTAG"
  print "testing if %s == ACT-CTAGCTAG" % aln2[1].seq       
  assert str(remove_dots(align1)[1].seq) == "ACT-CTAGCTAG"
  print "testing if %s == ACTGCAAGCTAG" % aln2[2].seq       
  assert str(remove_dots(align1)[2].seq) == "ACTGCAAGCTAG"


def test_remove_dots2():
  """test if sequence names are conserved"""
  aln2 = remove_dots(align1)
  print "testing if %s == Alpha" % aln2[0].id       
  assert aln2[0].id == "Alpha"
  print "testing if %s == Beta" % aln2[1].id       
  assert aln2[1].id == "Beta"
  print "testing if %s == Gamma" % aln2[2].id       
  assert aln2[2].id == "Gamma"

@raises(AssertionError)  
def test_remove_dots3():
  """test if error handling works correctly when dot in reference seq"""
  align1 = MultipleSeqAlignment([
             SeqRecord(Seq("A.TGCTAGCTAG", generic_dna), id="Alpha"),
             SeqRecord(Seq("ACT-CTAGC.AG", generic_dna), id="Beta"),
             SeqRecord(Seq("ACTGCTAGCTAG", generic_dna), id="Gamma"),
         ])
  remove_dots(align1)
  
def test_write_phylip():
  output = StringIO()
  write_phylip(align1, output)
  results = output.getvalue().split('\n')
  output.close()
  assert results[0] == '3 12'
  assert results[1] == 'Alpha     ACTGCTAGCTAG'
  assert results[2] == 'Beta      ACT-CTAGC.AG'
  assert results[3] == 'Gamma     .....A......'

def test_write_phylip2():
  """ if writer deals with long ids correctly"""
  align1[0].id = 'very_long_name'
  output = StringIO()
  write_phylip(align1, output)
  results = output.getvalue().split('\n')
  output.close()
  assert results[0] == '3 12'
  assert results[1] == 'very_long_name ACTGCTAGCTAG'
  assert results[2] == 'Beta           ACT-CTAGC.AG'
  assert results[3] == 'Gamma          .....A......'

def test_six_frame_translation():
  coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
  translations = six_frame_translation(coding_dna)
  print "test if %s == MAIVMGR*KGAR*" % str(translations[0])
  assert str(translations[0]) == "MAIVMGR*KGAR*"
  print "test if %s == LSGTLSAAHYNGH" % str(translations[1])
  assert str(translations[1]) == "LSGTLSAAHYNGH"
  print "test if %s == WPL*WAAERVPD" % str(translations[2])
  assert str(translations[2]) == "WPL*WAAERVPD"
  print "test if %s == YRAPFQRPITMA" % str(translations[3])
  assert str(translations[3]) == "YRAPFQRPITMA"
  print "test if %s == GHCNGPLKGCPI" % str(translations[4])
  assert str(translations[4]) == "GHCNGPLKGCPI"
  print "test if %s == IGHPFSGPLQWP" % str(translations[5])
  assert str(translations[5]) == "IGHPFSGPLQWP"

def test_find_CaaX():
  coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGATGTGCTGCCCGATAG", IUPAC.unambiguous_dna)
  assert len(find_CaaX(coding_dna)) == 0
  coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGGTGTGCTGCCCGATAG", IUPAC.unambiguous_dna)
  assert len(find_CaaX(coding_dna)) == 1
  
