#!/Users/HeathOBrien/anaconda/bin/python
import fileinput, getopt, sys, re, os
import pandas as pd

"""default filetype is transactions. graph not made by default"""

def main(argv):
    infilename = ''
    filetype = 'transactions'
    graph = 0
    balance_forward = ''
    usage = "HSBC.py -b balance_forward -i new_transactions [-t statement|transactions|csv] [ -g] >>~/Desktop/Desktop_files/transactions.txt"
    try:
        opts, args = getopt.getopt(argv,"hgi:b:t:",["help", "graph", "infile=", "filetype=", "balance="])
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print "Help:"
            print usage
            sys.exit()
        elif opt in ("-g", "--graph"):
            graph = 1
        elif opt in ("-i", "--ifile"):
            infilename = arg
        elif opt in ("-b", "--balance"):
            balance_forward = arg
        elif opt in ("-t", "--filetype"):
            filetype = arg

    try:
        balance_forward = float(balance_forward)
    except:
        sys.exit(usage)
    if filetype == 'transactions':
        ParseTransactions(infilename, balance_forward)
    elif filetype == 'csv':
        ParseCSV(infilename, balance_forward)
    else:
        ParseStatement(infilename, balance_forward)
    if graph == 1:
        #print "Rscript ~/Documents/R/PlotBalance.R"
        os.system("Rscript ~/Documents/R/PlotBalance.R")
        
def ParseStatement(infilename, balance):
    date_re = re.compile('^D?(\d\d) ((Jan)|(Feb)|(Mar)|(Apr)|(May)|(Jun)|(Jul)|(Aug)|(Sep)|(Oct)|(Nov)|(Dec))')
    type_re = re.compile('((ATM)|(BP)|(CR)|(DD)|(DR)|(SO)|(TFR)|(VIS))')
    amount_re = re.compile('(\d*\.\d\d)')
    date = ''
    venue = ''
    transaction_type = ''
    amount = 0
    linenum = 0
    for line in fileinput.input(infilename):
        linenum += 1
        line = line.strip()
        if not line:
            continue
        m = date_re.match(line)
        t = type_re.match(line)
        a = amount_re.match(line)
        if m:
            if date:
                sys.exit("line %s (%s) of transaction on %s formatted incorrectly" % (linenum, line, date))
            if m.group(2) == 'Dec':
                date = "%s-%s-14" % (m.group(1), m.group(2))
            else:
               date = "%s-%s-15" % (m.group(1), m.group(2))
        elif t:
            if transaction_type:
                sys.exit("line %s (%s) of transaction on %s formatted incorrectly" % (linenum, line, date))   
            transaction_type = t.group(1)
        elif a:
            if transaction_type or venue or date:  #skip amount matches after first (sometime there is a balance line)
                if not transaction_type or not venue:
                    sys.exit("lines above %s (%s) of transaction on %s formatted incorrectly" % (linenum, line, date))
                amount = a.group(1)
                if transaction_type == 'CR':
                    balance += float(amount)
                else:
                    balance -= float(amount)
                
                print '\t'.join((date, transaction_type, venue, amount, str(balance)))
                date = ''
                venue = ''
                transaction_type = ''
        else:
            venue += line.replace("'", "")

def ParseCSV(infilename, balance):
    transactions = pd.DataFrame.from_csv(infilename, sep=',', index_col=False)
    transactions = transactions.iloc[::-1]
    transactions['Balance'] = transactions.Amount.cumsum() + balance
    transactions['Type'] = ''
    transactions['Date'] = transactions['Date'].str.replace(' ', '-').str.replace('20', '')
    transactions = transactions[['Date', 'Type', 'Description', 'Amount', 'Balance']]
    print pd.DataFrame.to_csv(transactions, index=False, header=False, sep='\t')
    
def ParseTransactions(infilename, balance_forward):      
    first_line = re.compile('^D?(\d\d) ((Jan)|(Feb)|(Mar)|(Apr)|(May)|(Jun)|(Jul)|(Aug)|(Sep)|(Oct)|(Nov)|(Dec)) *\t *((ATM)|(BP)|(CR)|(DD)|(DR)|(SO)|(TFR)|(VIS)) *\t *(.*)')
    last_line = re.compile('^([\t ]+)(\d+\.\d\d)([\t ]+)(\d+\.\d\d)')
    date = ''
    venue = ''
    transaction_type = ''
    output = []
    linenum = 0
    for line in fileinput.input(infilename):
        linenum += 1
        line = line.rstrip()
        m = first_line.match(line)
        if m:
            if m.group(2) == 'Dec':
                date = "%s-%s-14" % (m.group(1), m.group(2))
            else:
                date = "%s-%s-15" % (m.group(1), m.group(2))
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
                venue = venue.replace('#', '') 
                output.append([date, transaction_type, venue, str(balance_change), str(balance)])
            elif '\t' in line:
                sys.exit("line %s (%s) of transaction on %s formatted incorrectly" % (linenum, line, date))
            else:
                venue = ' '.join((venue, line))

    reversed = output[::-1]
    output = []
    for line in reversed:
        balance_forward += float(line[-2])
        line[-1] = str(balance_forward)
        print '\t'.join(line)

if __name__ == "__main__":
   main(sys.argv[1:])

