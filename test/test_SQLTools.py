from SQLTools import combine_clusters
import MySQLdb as mdb
from nose.tools import *

def test_combine_clusters():
  con = mdb.connect('localhost', 'root', '', 'test')
  with con:
    cur = con.cursor()
    cur.execute("""DROP TABLE IF EXISTS `cluster_num`;""")
    
    cur.execute("""CREATE TABLE `cluster_num` (
                     `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
                     `seqid` char(25) DEFAULT NULL,
                     `cluster` int(11) unsigned NOT NULL,
                     PRIMARY KEY (`id`),
                     KEY `seqid_inx` (`seqid`),
                     KEY `cluster_inx` (`cluster`)
                   ) ENGINE=InnoDB DEFAULT CHARSET=latin1;"""
                )

    cur.execute("""INSERT INTO `cluster_num` (`id`, `seqid`, `cluster`)
                     VALUES
                	   (1,'Seq1',1),
                	   (2,'Seq2',1),
                	   (3,'Seq3',1),
                	   (4,'Seq1',2),
                	   (5,'Seq2',2),
                	   (6,'Seq4',2);"""
                )

  combine_clusters(2,1,'test')
  con = mdb.connect('localhost', 'root', '', 'test')
  with con:
    cur.execute("SELECT seqid, cluster FROM cluster_num")
    results = cur.fetchall()
    print results
    assert results == (('Seq1', 1), ('Seq2', 1), ('Seq3', 1),('Seq4', 1))