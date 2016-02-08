#!/usr/local/bin/python

import MySQLdb as mdb

def combine_clusters(cluster1, cluster2, db='test'):
  """replace the cluster number of cluster1 sequences with cluster2 (cluster1 sequences
     that are also in cluster2 are deleted"""
  con = mdb.connect('localhost', 'root', '', db)
  with con:
    cur = con.cursor()
    cur.execute("SELECT seqid FROM cluster_num WHERE cluster = %s", cluster1)
    for seq in cur.fetchall():
      cur.execute("SELECT seqid FROM cluster_num WHERE cluster = %s AND seqid = %s", (cluster2, seq[0]))
      if len(cur.fetchall()) > 0:
        #print "DELETE FROM cluster_num WHERE seqid = %s AND cluster = %s" % (seq[0], cluster1)
        cur.execute("DELETE FROM cluster_num WHERE seqid = %s AND cluster = %s", (seq[0], cluster1))
      else:
        cur.execute("UPDATE cluster_num SET cluster = %s WHERE seqid = %s", (cluster2, seq[0]))

def get_seqs(clusternum, db='Selaginella'):
  """"""
  seqs = []
  cur.execute("SELECT seqid, sequence FROM Sequences WHERE repseq = 'NR' AND clusternum = %s", clusternum)
  for (seqid, sequence) in cur.fetchall():
    seq_record = SeqRecord(Seq(sequence), id=seqid, description = '')
    if seqid[:4] in ('KRAU', 'MOEL', 'UNCc', 'WILD'):  #BLUELEAF sequence. Need to remove non-coding portions
      cur.execute("SELECT start, end, strand, gene_id FROM orfs WHERE seqid = %s", seqid)
      for (start, end, strand, gene_id) in cur.fetchall():
        strand = int(strand)        
        #Check that gene_id is properly formatted (with cluster 
        if gene_id.split('_')[0] not in seqid or gene_id.split('_')[1][0] != 'c': #gene_id not properly formated
          print "skipping %s (gene_id %s)" % (seqid, gene_id)
          continue
        if int(gene_id.split('_')[1][1:]) != clusternum: #gene_id properly formatted, but not in agreement with clusternum
          print "skipping %s (gene_id %s)" % (seqid, gene_id)
          continue        
        seq_record = seq_record[start-1:end]
        if strand == 1:
          name = '_'.join(map(str, [seqid, start, end]))
        elif strand == 0:
          name = '_'.join(map(str, [seqid, end, start]))
          seq_record = seq_record.reverse_complement()
          seq_record.description = ''
        else:
          sys.exit('strand %s not recognized' % strand)
        seq_record.id = name 
        #seq_record.description = ''
        seqs.append(seq_record)
    elif seqid[:3] in ('EFJ', 'Sel', 'ATH'):  #other non 1KP sequences:
      seqs.append(seq_record)
  cur.execute("SELECT sequences.seqid, sequences.sequence FROM cluster_num, sequences WHERE sequences.seqid = cluster_num.seqid AND cluster_num.cluster = %s", clusternum)
  for (seqid, sequence) in cur.fetchall():
    seq_record = SeqRecord(Seq(sequence), id=seqid, description = '')
    seqs.append(seq_record)
    
    
  return seqs
