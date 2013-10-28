#!/usr/local/bin/python

import sys, time, getopt
from os import path, system

def main(species, step, logfilename):
  step = int(step)
  #Make all of the file/folder names needed for pipeline
  raw_assembly = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Assemblies", species + "_Tr.fa")
  blast_all = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Blast", species + "_Tr_Selmo_all.bl")
  gtf_all = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "GTF", species + "_Tr_1kp.gtf")
  contig_folder = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "ContigClusters")
  consensus_file = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Assemblies", species + "_nr.fa")
  indexed_file = path.join(path.expanduser("~"), "Bioinformatics", "Index", species + "_Tr_cons.fa")
  index = path.join(path.expanduser("~"), "Bioinformatics", "Index", species + "_Tr_cons")
  merge_info = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Mappings", species + "_Tr_cons_unstranded_merge.txt")
  merge_dir = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Mappings", species + "_Tr_cons_unstranded_merge")
  consensus_blast = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Blast", species + "_Tr_cons_selmo_all.bl")
  consensus_blast_table = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Blast", species + "_Tr_cons_selmo_all_tbl.bl")

  if species == 'MOEL':
    lib_num = 3
  else:
    lib_num = 4
  
  #Create Log File to save commands
  log = open(logfilename, "a", 1)

  #Change sequence names
  if step <= 1:
    log.write("%s, Change sequence names:\n" % time.asctime(time.localtime()))
    log.write("perl -pi -e 's/>(comp\S+).*/>%s$1/g' %s\n" % ( species, raw_assembly ))
    system("perl -pi -e 's/>(comp\S+).*/>%s$1/g' %s\n" % ( species, raw_assembly ))

  #Blast all contigs agains database containing all Selmo 1kp sequences and the Ensembl annotation of S_moellendorffii 
  if step <= 2:
    log.write("\n%s, Blast all contigs against all Selaginella AA sequences:\n" % time.asctime(time.localtime()))
    log.write(" ".join(("blastx -query %s" % raw_assembly,
                "-db Selmo_all_aa", 
                "-out %s" % blast_all, 
                "-evalue 1e-20", 
                "-num_threads 8", 
                "-outfmt '6 qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore'\n"
                ))
             )
    system(" ".join(("blastx -query %s" % raw_assembly,
                "-db Selmo_all_aa", 
                "-out %s" % blast_all, 
                "-evalue 1e-20", 
                "-num_threads 8", 
                "-outfmt '6 qseqid qlen sacc slen pident length mismatch gapopen qstart qend qframe sstart send sframe evalue bitscore'"      
                ))
          )

  if step <= 4:
    log.write("\n%s, Add contig sequences to clusters:\n" % time.asctime(time.localtime()))
    log.write("Blast2OrthologGroup.py -b %s -s %s\n" % (blast_all, raw_assembly)) #Add sequences to clusters
    system("Blast2OrthologGroup.py -b %s -s %s" % (blast_all, raw_assembly)) #Add sequences to clusters

  #build trees for each ortholog group.
  if step <=5:
    log.write("\n%s, Build trees for each ortholog group:\n" % time.asctime(time.localtime()))
    log.write("ParallelBatch.py -f %s -g .fa -c 'RunPhyml.py * fast'\n" % contig_folder)
    system("ParallelBatch.py -f %s -g .fa -c 'RunPhyml.py * fast'" % contig_folder)










  #Make consensus sequences
  if step <=6:
    log.write("\n%s, Make consensus sequences:\n" % time.asctime(time.localtime()))
    log.write("rm names.p\n")
    system("rm names.p")
    log.write("ParallelBatch.py -p 1 -f %s -g .fa -c 'RunMakeConsensus.py * %s'\n" % (contig_folder, species) )
    system("ParallelBatch.py -p 1 -f %s -g .fa -c 'RunMakeConsensus.py * %s'" % (contig_folder, species) )
    #remove contig sequences from alignments
    log.write("ParallelBatch.py -p 8 -f %s -g .fa -c 'RemoveSeqs.py * %s'\n" % ( contig_folder, species) )
    system("ParallelBatch.py -p 8 -f %s -g .fa -c 'RemoveSeqs.py -s * -g %s'" % ( contig_folder, species) )
    
  #Index consensus file for read mapping
  if step <=7:
    log.write("\n%s, Index consensus file for read mapping:\n" % time.asctime(time.localtime()))
    log.write("cp %s %s\n" % (consensus_file, indexed_file))
    system("cp %s %s" % (consensus_file, indexed_file))
    log.write("bowtie2-build %s %s\n" % (indexed_file, index))
    system("bowtie2-build %s %s" % (indexed_file, index))

  #Run TopHat mapping and cufflinks on each library
  if step <=8:
    for x in range(lib_num):
      x = x + 1
  
      mapping_dir = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Mappings", species + str(x) + "_Tr_cons_unstranded")
      f_reads = ",".join((path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Filtered", species + str(x) + "_1.filtered.fastq"), 
                          path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Filtered", species + str(x) + "b_1.filtered.fastq")
                        ))
      r_reads = ",".join((path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Filtered", species + str(x) + "_2.filtered.fastq"), 
                          path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Filtered", species + str(x) + "b_2.filtered.fastq")
                        ))
      accepted_hits = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Mappings", species + "_Tr_cons_unstranded", "accepted_hits.bam")
      bam_file = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Mappings", species + "_Tr_cons_unstranded", species + str(x) + "Tr_cons.bam")
      cufflinks_folder = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Mappings", species + "_Tr_cons_unstranded", "Cufflinks")
      cufflinks_gtf = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Mappings", species + "_Tr_cons_unstranded", "Cufflinks", "transcripts.gtf")

      log.write("\n%s, Map %s%s reads consensus sequences:\n" % (time.asctime(time.localtime()), species, str(x)))
      log.write(" ".join(("tophat2 --library-type fr-unstranded",
                  "--num-threads 8",
                  "--max-multihits 1",
                  "--output-dir %s" % mapping_dir,
                  index,
                  f_reads,
                  r_reads,
                  "\n"
                  ))
               )
      system(" ".join(("tophat2 --library-type fr-unstranded",
                  "--num-threads 8",
                  "--max-multihits 1",
                  "--output-dir %s" % mapping_dir,
                  index,
                  f_reads,
                  r_reads
                  ))
            )
      log.write("\n%s, Index %s%s read mapping:\n" % (time.asctime(time.localtime()), species, x))
      log.write("mv %s %s\n" % (mapping_dir + "/accepted_hits.bam", bam_file))
      system("mv %s %s" % (mapping_dir + "/accepted_hits.bam", bam_file))
  
      log.write("samtools index %s\n" % bam_file)
      system("samtools index %s" % bam_file)

      log.write("\n%s, Run cufflinks of %s%s read mapping:\n" % (time.asctime(time.localtime()), species, x))  
      log.write("cufflinks -o %s %s\n" % ( cufflinks_folder, bam_file))
      system("cufflinks -o %s %s" % ( cufflinks_folder, bam_file))
   
      mergeinfo_fh = open(merge_info, "a")
      mergeinfo_fh.write("%s\n" % cufflinks_gtf)
      mergeinfo_fh.close()
  
  #Run cuffmerge on transcripts
  if step <=9:
    log.write("\n%s, Run cuffmerge on cufflinks transcripts:\n" % time.asctime(time.localtime()))
    log.write("cuffmerge -o %s -s %s %s\n" % (merge_dir, consensus_file, merge_info)) 
    system("cuffmerge -o %s -s %s %s" % (merge_dir, consensus_file, merge_info))

  #This is as far as I go so long as I need to use the online version of OrfPredictor
  log.close()


if __name__ == "__main__":
  step = 7
  species = ''
  logfilename = ''
  try:
    opts, args = getopt.getopt(argv,"hn:s:l:",["name=", "step=","log="])
  except getopt.GetoptError:
    print 'SelaginellaPipeline.py -n <name> [ -s <step> -l <logfile> ]'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print 'SelaginellaPipeline.py -n <name> [ -s <step> -l <logfile> ]'
      sys.exit()
    elif opt in ("-n", "--name"):
      species = arg
    elif opt in ("-l", "--log"):
      logfilename = arg
    elif opt in ("-s", "--step"):
      step = arg
  if not logfilename:
    logfilename = path.join(path.expanduser("~"), "Bioinformatics", "Selaginella", "Logs", species + "_log.txt")    
  if species:
    main(species, step, logfilename)  
  else:  
    for species in ('UNC', 'MOEL', 'WILD', 'KRAUS'):
      main(species, step, logfilename)
  exit