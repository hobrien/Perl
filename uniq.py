#!/usr/local/bin/python


import argparse, sys

parser = argparse.ArgumentParser(description='like unix uniq, but allows specification of a specific field and need not be sorted')
parser.add_argument('files', metavar='files', type=argparse.FileType('r'), nargs='?',
                   default=sys.stdin, help='list of files to process')
parser.add_argument('-k, --column', dest='column', default=1, type=int,
                   help='column used for comparison of unique')
parser.add_argument('-d, --delimitor', dest='sep', default='', type=str,
                   help='column used for comparison of unique')

args = parser.parse_args()
if args.sep == '\\t' or args.sep == 't':
  args.sep = r'\t'.decode('string-escape')

previous = {}
for line in args.files:
  line = line.strip()
  if args.sep:
    current = line.split(args.sep)[args.column - 1]
  else:
    current = line.split()[args.column - 1]
  if current not in previous:
    print line
    previous[current] = 1
