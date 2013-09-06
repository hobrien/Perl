#!/opt/local/bin/python

import sys, os, time
from subprocess import Popen, list2cmdline
import glob
import sys, getopt
from os import path, system

def main(argv):
  command = ''
  searchterm = '' #limit to files that match search term
  folder = ''
  try:
      opts, args = getopt.getopt(argv,"hc:g:f:",["command=","grep=", "folder="])
  except getopt.GetoptError:
    print 'Type BatchCommand.py -h for options'
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-h":
       print 'BatchCommand.py -c <command> -g <searchterm>'
       sys.exit()
    elif opt in ("-c", "--command"):
       command = arg
    elif opt in ("-g", "--grep"):
       searchterm = arg
    elif opt in ("-f", "--folder"):
       folder = arg
  commands = []
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

    return num - 2

def exec_commands(cmds):
    ''' Exec commands in parallel in multiple process 
    (as much as we have CPU)
    '''
    if not cmds: return # empty list

    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)

    max_task = cpu_count()
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