import sys
import re

'''
In sam2tsv output
 read position is 0-based
 reference posiiton is 1-based
input: sam2tsv file
#READ_NAME      FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
'''

if len(sys.argv) != 2:
    print "\n usage: python single_site_var.stats.py tsv_file \n" #ref   posiiton\n"
    exit(0)

sam2tsv = sys.argv[1]
#ref = sys.argv[2]
#position = sys.argv[3]

rdnames = []
qualities = {}
dels = []
mis = []
ins = []
match = []
ref_pos = {}

with open (sam2tsv,'r') as fh:
    for line in fh:
        if re.match ('\s+',line):
            continue
        if re.match('#',line):
            continue
        if re.match (':',line):
            continue

        ary = line.strip().split()

        if len (ary) != 9:
            continue
        if ary[6].startswith('-') : #ref position
            continue
#READ_NAME      FLAG    CHROM   READ_POS  BASE  QUAL   REF_POS REF   OP
        if ary[6] != '.':
            ary[6] = int (ary[6])   # ref pos is 1-based
        if ary[3] != '.':
            ary[3] = int (ary[3]) + 1  # read pos is 0-based
        k = (ary[2], ary[6], ary[7], ary[0], ary[3],ary[4])
        qualities[k] = ord(ary[5]) - 33
        if ary[-1].upper() == 'M' and ary[4] == ary[7]: # or re.match ('I', ary[-1].upper()):
            match.append(k)
            rdnames.append(k)
        if ary[-1].upper() == 'M' and ary[4] != ary[7]: # or re.match ('I', ary[-1].upper()):
            mis.append(k)
            rdnames.append(k)
        elif ary[-1].upper() == 'I':
            ins.append(k)
            rdnames.append(k)
        elif ary[-1].upper() == 'D':
            dels.append(k)
            rdnames.append(k)
        else:
            continue

print ",".join (["#REF",'REF_POS','REF_BASE','READ_NAME','READ_POSITION','READ_BASE','BASE_QUALITY','MISMATCH','INSERTION','DELETION'])

for k in rdnames:
    #print ref_pos[k], k, qualities[k],
    new_k = ','.join (map (str, list (k)))
    q = str (qualities[k])
    if k in mis:
        m =  '1'
    else:
        m = '0'
    if k in ins:
        i = '1'
    else:
        i = '0'
    if k in dels:
        d = '1'
    else:
        d = '0'
    print ','.join ([new_k,q,m,i,d])

