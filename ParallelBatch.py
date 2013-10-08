#!/usr/local/bin/python

"""take all files in specified folder that match the search criteria and replace asterisks
in the command with the filenames and execute it on the specified number of processors 
(default = 6)
"""
import sys, os, time
from subprocess import Popen, list2cmdline
import glob
import sys, getopt
from os import path, system


def main(argv):
  command = ''
  searchterm = '' #limit to files that match search term
  folder = ''
  max_task = 6
  try:
      opts, args = getopt.getopt(argv,"hc:g:f:p:",["command=","grep=", "folder=", "processors="])
  except getopt.GetoptError:
    print 'Type ParallelBatch.py -h for options'
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print 'ParallelBatch.py -c <command> -g <searchterm> -f <folder> -p <num_procesors>'
       sys.exit()
    elif opt in ("-c", "--command"):
       command = arg
    elif opt in ("-g", "--grep"):
       searchterm = arg
    elif opt in ("-f", "--folder"):
       folder = arg
    elif opt in ("-p", "--processor"):
       max_task = arg

  commands = [max_task]
  for file in glob.glob(path.join(folder, '*')):
    if searchterm == "" or searchterm in file: 
      commands.append(command.replace('*',file).split(" "))
  
  exec_commands(commands)

def cpu_count():
    ''' Returns the number of CPUs in the system
    '''
    num = 1
    if sys.platform == 'win32':
        try:
            num = int(os.environ['NUMBER_OF_PROCESSORS'])
        except (ValueError, KeyError):
            pass
    elif sys.platform == 'darwin':
        try:
            num = int(os.popen('sysctl -n hw.ncpu').read())
        except ValueError:
            pass
    else:
        try:
            num = os.sysconf('SC_NPROCESSORS_ONLN')
        except (ValueError, OSError, AttributeError):
            pass

    return num

def exec_commands(cmds):
    ''' Exec commands in parallel in multiple process 
    (as much as we have CPU)
    '''
    max_task = int(cmds.pop(0))
    if not cmds: return # empty list

    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)

    processor_num = cpu_count()
    if max_task > processor_num:
      sys.exit("maximum number of processors is %i" % processor_num)
    processes = []
    while True:
        while cmds and len(processes) < max_task:
            task = cmds.pop()
            print list2cmdline(task)
            processes.append(Popen(task))

        for p in processes:
            if done(p):
                if success(p):
                    processes.remove(p)
                else:
                    fail()

        if not processes and not cmds:
            break
        else:
            time.sleep(0.05)


  
if __name__ == "__main__":
   main(sys.argv[1:])