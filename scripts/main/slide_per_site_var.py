import sys
import re
from collections import OrderedDict

mem = {}
contents = OrderedDict()

win_size = int (sys.argv[2])
dist = int (win_size) // 2 + 1

with open (sys.argv[1],'r') as fh:
    #Ref,pos,base,cov,q_mean,q_median,q_std,mis,ins,del
    for line in fh:
        if line.startswith ('#'):
            continue
        if re.match ('\s+',line):
            continue
        ary = line.strip().split(',')
        ref, pos = ary[0], ary[1]
        contents[ref+':'+pos] = line.strip()

header = '#Kmer,window,Relative_pos,Ref,Ref_Pos,base,cov,q_mean,q_median,q_std,mis,ins,del'
print header
for k in contents.keys():
    ref, pos = k.split(':')
    try:
	pos = int(pos)
    except:
	sys.stderr.write(k)
	continue

    upper = ''
    down = ''
    POS = []
    LINES = []
    for i in list (reversed (range (1,dist))):
        POS.append (str(pos-i))
        kk = ref+':'+str(pos-i)
        if kk in contents:
            base = contents[kk].split(',')[2]
            upper += base
            LINES.append ('-'+str(i)+','+ contents[kk])
        else:
            upper += 'N'
            LINES.append ('-'+str(i)+','+'Null')
    LINES.append ('+0'+','+contents[k])
    POS.append (str(pos))
    for j in range (1,dist):
        POS.append (str(pos + j))
        kk = ref+':'+str(pos + j)
        kk = ref+':'+str(pos+j)
        if kk in contents:
            base = contents[kk].split(',')[2]
            down += base
            LINES.append ('+'+str(j)+','+ contents[kk])
        else:
            down += 'N'
            LINES.append ('+'+str(j)+','+'None')
    positions = ':'.join (POS)
    for l in LINES:
        print upper+contents[k].split(',')[2] +down+',' +positions+','+ l
