import fileinput, sys, re

balance_forward = float(sys.argv[1])
first_line = re.compile('^D?(\d\d) ((Jan)|(Feb)|(Mar)|(Apr)|(May)|(Jun)|(Jul)|(Aug)|(Sep)|(Oct)|(Nov)|(Dec)) *\t *((ATM)|(BP)|(CR)|(DD)|(DR)|(SO)|(TFR)|(VIS)) *\t *(.*)')
last_line = re.compile('^([\t ]+)(\d+\.\d\d)([\t ]+)(\d+\.\d\d)')
date = ''
venue = ''
transaction_type = ''
output = []
for line in fileinput.input(sys.argv[2:]):
  line = line.rstrip()
  m = first_line.match(line)
  if m:
    date = "%s-%s-14" % (m.group(1), m.group(2))
    transaction_type = m.group(15)
    venue = m.group(24)
  else:
    m = last_line.match(line)
    if m:
      if m.group(1).count('\t') == 2 and m.group(3).count('\t') == 1:
        balance_change = float(m.group(2))
      elif m.group(1).count('\t') == 1 and m.group(3).count('\t') == 2:
        balance_change = 0.00 - float(m.group(2))
      else:
        sys.exit("line %s of transaction on %s formatted incorectly" % (line, date))
      balance = m.group(4)
      venue = venue.replace('UNIV OF B EXP AC', 'BRISTOL UNI')
      venue = venue.replace('PAYROLL A/C', 'BATH UNI')      
      venue = venue.replace('PROPERTY CONCEPT OBRIEN DIEZMANN', 'PROPERTY CONCEPT') 
      venue = venue.replace("'", '') 
      venue = venue.replace('"', '') 
      output.append([date, transaction_type, venue, str(balance_change), str(balance)])
    elif '\t' in line:
      sys.exit("line %s of transaction on %s formatted incorrectly" % (line, date))
    else:
      venue = ' '.join((venue, line))

reversed = output[::-1]
output = []
for line in reversed:
  balance_forward += float(line[-2])
  line[-1] = str(balance_forward)
  print '\t'.join(line)


