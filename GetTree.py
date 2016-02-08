#!/usr/local/bin/python

import sys, getopt
import MySQLdb as mdb
from os import path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import Queue, threading, time


#Append trimal results to Nexus file as an exclusion set
def add_exset(file, trimal_res):
  import re
  from itertools import groupby, count
  incl = trimal_res.split(',')
  incl = map(int, incl)
  aln = open(file, 'r')
  for line in aln.readlines():
    match = re.search( r"nchar *= *(\d+);", line, re.I)
    if match:
      length = int(match.group(1))
      break
  print length
  aln.close()
  excl = []
  for x in range(length):
    if not x in incl:
      excl.append(x+1)

  aln = open(file, 'a')
  aln.write("\nBEGIN ASSUMPTIONS;\n\t EXSET * TRIMAL  =  ")
  aln.write(' '.join(as_range(g) for _, g in groupby(excl, key=lambda n, c=count(): n-next(c))))  
  aln.write(";\n\nEND;\n")
  aln.write("BEGIN CODONS;\n\tCODONPOSSET * UNTITLED  =  1: 1 - %s\\3, 2: 2 - %s\\3, 3: 3 - %s\\3;\n\nEND;\n" % (length, length, length))
  aln.close

#convert trimal intervals to ranges   
def as_range(iterable):
  l = list(iterable)
  if len(l) > 1:
    return '{0}-{1}'.format(l[0], l[-1])
  else:
    return '{0}'.format(l[0])

class myThread (threading.Thread):
    def __init__(self, threadID, name, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
    def run(self):
        process_data(self.name, self.q)

def process_data(threadName, q):
    while not exitFlag:
        queueLock.acquire()
        if not workQueue.empty():
            (dirname, clusternum, num_seqs) = q.get()
            queueLock.release()
            nexus_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.nex')
            phy_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.phy')
            tree_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.nwk')
            pdf_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.pdf')
            if num_seqs < 4:
              print "cluster %s contains %s seqs. At least 4 required for tree bulding" % (clusternum, num_seqs)
            else:  
              print "ConvertAln.py -f phylip -i %s" % nexus_file
              subprocess.call(["ConvertAln.py -f phylip -i " + nexus_file], shell=True)
              if num_seqs < 50:
                print "phyml  --quiet --no_memory_check -i %s" % phy_file
                subprocess.call(["phyml --quiet --no_memory_check -i " + phy_file], shell=True)
              else:
                print "phyml  --quiet --no_memory_check -o n -b 0 -i %s" % phy_file
                subprocess.call(["phyml --quiet --no_memory_check -o n -b 0 -i " + phy_file], shell=True)      
              subprocess.call(["rm", phy_file])
              subprocess.call(["mv", phy_file + '_phyml_tree.txt', tree_file])
              subprocess.call(["rm", phy_file + '_phyml_stats.txt'])
              print "ColourTree.py -t %s -o %s" % (tree_file, pdf_file)
              subprocess.call(["ColourTree.py", "-t", tree_file, "-o", pdf_file])

        else:
            queueLock.release()
        time.sleep(1)
  

if __name__ == "__main__":
  exitFlag = 0
  dirname = ''
  clusternum = ''
  first = 1
  last = 53127
  usage = "GetTree.py -d <dirname> -c <cluster> | ( -f <first> -l <last> )"
  try:
     opts, args = getopt.getopt(sys.argv[1:],"hd:f:l:c:",["cluster=", "dir=", "first=", "last="])
  except getopt.GetoptError, e:
     print e
     print usage
     sys.exit(2)
  for opt, arg in opts:
     if opt == '-h':
        print usage
        sys.exit()
     elif opt in ("-c", "--cluser"):
        clusternum = int(arg)
     elif opt in ("-f", "--first"):
        first = int(arg)
     elif opt in ("-l", "--last"):
        last = int(arg)
     elif opt in ("-d", "--dir"):
        dirname = arg

  con = mdb.connect('localhost', 'root', '', 'test')
  with con:
    cur = con.cursor()
    try:
      first = clusternum - 1
      last = clusternum + 1
    except TypeError:
      first -= 1
      last += 1  
    cur.execute("SELECT cluster FROM OrthoMCL WHERE cluster > %s AND cluster < %s GROUP BY cluster", (first, last))
    clusters = []
    for row in cur.fetchall():
      clusters.append(row[0])
    num_seqs = []
    print "Writing Sequences to file"
    for clusternum in clusters:
      seqs = []
      cur.execute("SELECT Sequences.seqid, Sequences.sequence FROM OrthoMCL, Sequences WHERE Sequences.seqid = OrthoMCL.seqid AND OrthoMCL.cluster = %s", clusternum)
      for (seqid, sequence) in cur.fetchall():
        seq_record = SeqRecord(Seq(sequence), id=seqid, description = '')
        seqs.append(seq_record)
      seq_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.fa')
      SeqIO.write(seqs, seq_file, "fasta")
      num_seqs.append(len(seqs))
    
    print "making alignments"
    index = -1
    for clusternum in clusters:
      index += 1
      if num_seqs[index] < 2:
        print "cluster %s contains %s seqs. At least 2 required for alignment" % (clusternum, num_seqs[index])
        continue
      seq_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.fa')
      translatorx_file = path.join(dirname, 'translatorx_res.nt_ali.fasta')
      aln_file = path.join(dirname, 'Cluster_' + str(clusternum) + '_aln.fa')
      nexus_file = path.join(dirname, 'Cluster_' + str(clusternum) + '.nex')
      print "translatorx_vLocal.pl -i %s -p F" % seq_file
      subprocess.call(["translatorx_vLocal.pl", "-i", seq_file, "-p", "F"])
      subprocess.call(["ConvertAln.py", "-i", translatorx_file, "-o", nexus_file, "-f", "nexus"])
      proc = subprocess.Popen(["trimal -gappyout -in %s -out %s -colnumbering" % (translatorx_file, aln_file)], stdout=subprocess.PIPE, shell=True)
      (trimal_res, err) = proc.communicate()
      add_exset(nexus_file, trimal_res)
      subprocess.call(["rm " + "translatorx_*"], shell=True)
      #subprocess.call(["rm", aln_file])

    print "Making Trees"
    threadList = ["Thread-1", "Thread-2", "Thread-3", "Thread-4", "Thread-5", "Thread-6"]
    queueLock = threading.Lock()
    workQueue = Queue.Queue(0)
    threads = []
    threadID = 1

    # Create new threads
    for tName in threadList:
      thread = myThread(threadID, tName, workQueue)
      thread.start()
      threads.append(thread)
      threadID += 1

    # Fill the queue
    queueLock.acquire()
    index = -1
    for clusternum in clusters:
      index +=1
      workQueue.put((dirname, clusternum, num_seqs[index]))
    queueLock.release()

    # Wait for queue to empty
    while not workQueue.empty():
      pass
    
    # Notify threads it's time to exit
    exitFlag = 1

    # Wait for all threads to complete
    for t in threads:
      t.join()


