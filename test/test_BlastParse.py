
from types import *
import os, warnings, sys
from nose.tools import *
from nose import with_setup
from ParseBlast import parse_args, parse_blast_stats, top_query, get_seq, combine_hits, parse_blast

sys.path.append(os.path.join(os.path.expanduser('~'), 'Perl'))

"""This is to be used with nosetests: 'nosetests ~/Perl' """
test_dir = os.path.join(os.path.expanduser('~'), 'Perl', 'test')


  
def test_parse_blast_stats():
  column_names = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" 
  row = "jgi|Copci1|227	901175	61.77	361	89	5	1309	1623	68805	67732	0.0	  413"
  results = parse_blast_stats(column_names, row.split())
  assert len(results.keys()) == 13
  assert results['qseqid'] == 'jgi|Copci1|227'
  assert results['sseqid'] == '901175'
  assert type(results['sseqid']) is StringType
  assert results['pident'] == 61.77
  assert type(results['pident']) is FloatType
  assert results['length'] == 361
  assert type(results['length']) is IntType
  assert results['mismatch'] == 89
  assert type(results['mismatch']) is IntType
  assert results['gapopen'] == 5
  assert type(results['gapopen']) is IntType
  assert results['qstart'] == 1309
  assert type(results['qstart']) is IntType
  assert results['qend'] == 1623
  assert type(results['qend']) is IntType
  assert results['sstart'] == 68805
  assert type(results['sstart']) is IntType
  assert results['send'] == 67732
  assert type(results['send']) is IntType
  assert results['evalue'] == 0.0
  assert type(results['evalue']) is FloatType
  assert results['bitscore'] == 413
  assert type(results['bitscore']) is FloatType
  assert results['strand'] == -1
  assert type(results['strand']) is IntType

def test_parse_args():
  args = parse_args((os.path.join(test_dir, 'test.tblastn'), os.path.join(test_dir, 'test.fa')))
  assert args.blastfilename ==  os.path.join(test_dir, 'test.tblastn')
  assert args.seqfilename ==  os.path.join(test_dir, 'test.fa')
  assert args.column_names == 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
  assert args.evalue == 10
  assert args.index_filename == os.path.join(test_dir, 'test.inx')
  assert args.program == 'tblastn'
  assert args.subject_name == 'test'
  assert args.translate is False
  args = parse_args(('--translate', os.path.join(test_dir, 'test.tblastn'), os.path.join(test_dir, 'test.fa')))
  assert args.translate
  args = parse_args(('--translate', '--evalue', '1e-20', os.path.join(test_dir, 'test.tblastn'), os.path.join(test_dir, 'test.fa')))
  assert args.evalue == 1e-20

def test_top_query():  #I can add to this to test the evalue filtering
  args = parse_args((os.path.join(test_dir, 'test.tblastn'), os.path.join(test_dir, 'test.fa')))
  results = top_query(args)
  assert len(results) == 11
  assert len(results[0]) == 12
  assert results[0][0]['pident'] == 61.77
  
def test_get_seqs():
  seqname = '1'
  args = parse_args((os.path.join(test_dir, 'test.tblastn'), os.path.join(test_dir, 'test.fa')))
  seq = get_seq(args, seqname)
  print str(seq)
  assert str(seq) == 'CCCACACCGCGCAAAATTCTATACCAAAATCGACCT'    
  args = parse_args(('--translate', os.path.join(test_dir, 'test.tblastn'), os.path.join(test_dir, 'test.fa')))
  seq = get_seq(args, seqname, start=2)
  print str(seq)
  assert str(seq) == 'PHRAKFYTKID'    

def test_combine_hits():
  """this will ensure that all pairs of coordinates are in the correct order (lowest to
  highest for + strand hits; highest to lowest for - strand hits"""
  
  args = parse_args((os.path.join(test_dir, 'test.tblastn'), os.path.join(test_dir, 'test.fa')))
  results = top_query(args)
  for result in results:
    coords =  combine_hits(result)
    print coords
    if coords[0][0] < coords[0][1]:
      for x in range(len(coords)-1):
        assert coords[x][0] > coords[x+1][0]
    else:
      for x in range(len(coords)-1):
        assert coords[x][0] < coords[x+1][0]
    
def test_combine_hits_overlap():
  """If the overlap removal is working correctly, the combined length of the coordinates
  should be 3 * the total length of the query (6 bp less than combined length of subjects"""
  
  results = [{'qstart': 1, 'qend': 5, 'sstart': 1, 'send': 15, 'strand': 1},
             {'qstart': 4, 'qend': 8, 'sstart': 21, 'send': 35, 'strand': 1}]
  coords =  combine_hits(results)
  print coords
  assert coords[0][1]-coords[0][0]+1 + coords[1][1]-coords[1][0]+1 == (results[1]['qend'] - results[0]['qstart'] + 1) * 3
    
def test_parse_blast():
  """Putting it all together"""
  args = parse_args(('--translate', os.path.join(test_dir, 'test.tblastn'), os.path.join(test_dir, 'test.fa')))
  seqs = parse_blast(args)
  assert len(seqs) == 11
  print seqs[0].seq[:30]
  assert seqs[0].seq.find("RSQTPNSSVADDTRLHKRSMTISKGHSVSVVLISAALETIAASKEARRSTVIRESTQK") == 0
  print seqs[0].seq.count('*')
  assert seqs[0].seq.count('*') == 13
