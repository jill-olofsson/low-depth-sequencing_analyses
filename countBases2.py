import sys
import re

inFile = open(sys.argv[1],'r')

print 'chr\tpos\tref\tcov\tA\tC\tG\tT\tdel\tins\tinserted\tambiguous'

for line in inFile:
    data = line.strip().split('\t')
    chrm = data[0]
    pos = data[1]
    ref = data[2].upper()
    cov = data[3]
    possibleSNP = ""
    
    types = {'A':0,'C':0,'G':0,'T':0,'-':0,'+':[],'X':[]}
    
    if ( int(cov) > 0 ):
        bases = data[4].upper()

        i = 0
        while i < len(bases):
            base = bases[i]

            if base == '^' or base == '$':
                i += 1
            elif base == '-':
                i += 1
            elif base == '*':
                types['-'] += 1
            elif base == '+':
                match = re.search("(\d*).*",str(bases[i+1:]))
                if match:
                    addLen = len(match.group(1))
                    addNum = int(match.group(1))
                else:
                    print "Match error!"
                    exit(2)
                    
                i += addLen
                #addNum = int(bases[i]) # we hope there are at most 9 insertions
                addSeq = '' # insertion sequence
                for a in range(addNum):
                        i += 1
                        addSeq += bases[i]

                types['+'].append(addSeq)
            elif base == '.' or base == ',':
                types[ref] += 1
            else:
                if types.has_key(base):
                        types[base] += 1
                else:
                        types['X'].append(base)

            i += 1
    
    adds = '.'
    if len(types['+']) > 0:
        adds = ','.join(types['+'])
    
    amb = '.'
    if len(types['X']) > 0:
        amb = ','.join(types['X'])
    
    if max([types[k] for k in ['A','C','G','T']]) != types[ref]:
        possibleSNP = "*"
    
    out = [chrm,pos,''.join([ref,possibleSNP]),cov,types['A'],types['C'],types['G'],types['T'],types['-'],len(types['+']),adds,amb]
    print '\t'.join([str(x) for x in out])
