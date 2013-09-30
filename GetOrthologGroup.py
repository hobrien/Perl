#!/opt/local/bin/python
"""Use GTF file to add sequences to ortholog groups
"""
import sys, getopt
from os import listdir, path, system
from Bio import SeqIO
from BCBio import GFF
from StringIO import StringIO

def main(argv):
  gtf_filename = ''
  seqfilename = ''
  feature = 'cds'
  seqtype = ''
  try:
      opts, args = getopt.getopt(argv,"hg:s:f:t:",["seqfile=","blastfile=", "feature=", "seqtype="])
  except getopt.GetoptError:
    print 'Type GetOrthologGroups.py -h for options'
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print 'GetOrthologGroups.py -g <gtf_file> -s <seqfile> -f <feature> -t <seqtype>'
       sys.exit()
    elif opt in ("-g", "--gtffile"):
       gtf_filename = arg
    elif opt in ("-s", "--seqfile"):
       seqfilename = arg
    elif opt in ("-f", "--feature"):
       feature = arg
    elif opt in ("-t", "--seqtype"):
       seqtype = arg
  if seqtype != "contigs" or seqtype == "consensus":
    sys.exit("seqtype must be either 'contigs' or 'consensus'")     
  seqfilehandle = open(seqfilename)
  seq_dict = SeqIO.to_dict(SeqIO.parse(seqfilehandle, "fasta"))
  seqfilehandle.close()
  gtf_filehandle = open(gtf_filename)
  for SeqRec in GFF.parse(gtf_filehandle, base_dict=seq_dict): 
    if not SeqRec.features:
      continue                                                  #Skip sequences that are not in the GFF
    #cluster_num = SeqRec.features[0]
    cluster_num = SeqRec.features[0].id.split('_')[1]
    if not cluster_num[0] =='c':                               #Skip sequences that match reference sequences that are not part of a cluster
      continue
    cluster_num = cluster_num[1:]
    if seqtype == "contigs":
      cluster_filename = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "ContigClusters", "Cluster_" + cluster_num + ".fa")
    elif seqtype == "consensus":
      cluster_filename = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "ConsensusCLusters", "Cluster_" + cluster_num + ".fa")     
    if feature == 'cds':
      subseq = SeqRec.features[0].sub_features[1].extract(SeqRec)
    elif feature == 'contigs':
      if SeqRec.features[0].sub_features[1].location.strand == -1:
        subseq = SeqRec.seq.reverse_complement()
      else:
        subseq = SeqRec.seq 
    else:
      sys.exit("feature %s not recognized" % feature)
      
    if SeqRec.features[0].sub_features[1].location.strand == -1:
      subseq.id = SeqRec.id + '_rc'
      subseq.description = subseq.id
    cluster_file = open(cluster_filename, "a")
    cluster_file.write(subseq.format("fasta"))
    cluster_file.close()
  gtf_filehandle.close()
  
if __name__ == "__main__":
   main(sys.argv[1:])


