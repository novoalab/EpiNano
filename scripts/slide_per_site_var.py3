#!/usr/bin/env python3
import sys,os
import re
from collections import OrderedDict
from collections import defaultdict
import gzip

if len (sys.argv) < 2:
    print ("USAGE:\n  python " + sys.argv[0] + ' per_site_var ' )
    exit(0)


mem = {}
contents = OrderedDict()

win_size = 5
dist = int (win_size) // 2 + 1

per_site_var = sys.argv[1]

if per_site_var.endswith ('.gz'):
    fh = gzip.open (per_site_var,'r')
else:
    fh = open (per_site_var,'r')

#with open (sys.argv[1],'r') as fh:
    #Ref,pos,base,cov,q_mean,q_median,q_std,mis,ins,del
for line in fh:
    if line.startswith ('#'):
        continue
    if re.match ('\s+',line):
        continue
    ary = line.strip().split(',')
    ref, pos = ary[0], ary[1]
    contents[ref+':'+pos] = line.strip()
prefix = ".".join (sys.argv[1].split('.')[:-1])
out1 = prefix +'.per_site_var.5mer.tmp'
outh1= open (out1,'w')

header = '#Kmer,window,Relative_pos,Ref,Ref_Pos,base,cov,q_mean,q_median,q_std,mis,ins,del'
outh1.write (header+'\n')

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
        outh1.write (upper+contents[k].split(',')[2] +down+',' +positions+','+ l+'\n')
outh1.close()
############################################################################################
###    sum up slided per site variants from multiple lines into single lines ###############
############################################################################################
'''
input per site
#Kmer,window,Relative_pos,Ref,Ref_Pos,base,cov,q_mean,q_median,q_std,mis,ins,del
NNGTC,43:44:45:46:47,-2,Null
NNGTC,43:44:45:46:47,-1,Null
NNGTC,43:44:45:46:47,+0,cc6m_2244_T7_ecorv,45,G,1,5.000,5.000,0.000,0.000,0.000,0.000
NNGTC,43:44:45:46:47,+1,cc6m_2244_T7_ecorv,46,T,3,7.333,8.000,1.700,0.000,0.000,0.000
NNGTC,43:44:45:46:47,+2,cc6m_2244_T7_ecorv,47,C,3,5.333,6.000,1.700,0.000,0.000,0.000
NGTCG,44:45:46:47:48,-2,Null
NGTCG,44:45:46:47:48,-1,cc6m_2244_T7_ecorv,45,G,1,5.000,5.000,0.000,0.000,0.000,0.000
NGTCG,44:45:46:47:48,+0,cc6m_2244_T7_ecorv,46,T,3,7.333,8.000,1.700,0.000,0.000,0.000
NGTCG,44:45:46:47:48,+1,cc6m_2244_T7_ecorv,47,C,3,5.333,6.000,1.700,0.000,0.000,0.000
NGTCG,44:45:46:47:48,+2,cc6m_2244_T7_ecorv,48,G,3,10.000,9.000,1.414,0.000,0.000,0.000
GTCGA,45:46:47:48:49,-2,cc6m_2244_T7_ecorv,45,G,1,5.000,5.000,0.000,0.000,0.000,0.000
'''

mem_window = defaultdict (defaultdict(list).copy)
k_pool = []
f =  open (out1,'r')
for l in f:
    if l.startswith ('N'):
        continue
    elif l.startswith ('#'):
        continue
    ary = l.strip().split(',')
    if ary[0].endswith ('N'):
        continue
    ks = ary[0:4]
    del ks[2]
    k = ",".join (ks)
    k_pool.append (k)
    #print l.strip()
    if (len(ary) < 13):
        continue
    try:
        mem_window[k]['q'].append(ary[7])
    except:
        sys.stderr.write ("problematic line for q: " + l.strip())
    try:
        mem_window[k]['m'].append(ary[10])
    except:
        sys.stderr.write ("problematic line for m" + l.strip())
    try:
        mem_window[k]['i'].append(ary[11])
    except:
        sys.stderr.write ("problematic line for i:" + l.strip())
    try:
        mem_window[k]['d'].append(ary[12])
    except:
        sys.stderr.write ("problematic line for d" + l.strip())
    try:
        mem_window[k]['cov'].append(ary[6])
    except:
        sys.stderr.write ('problematic line for cov:' + l.strip())
f.close()
out2 = prefix + '.per_site_var.5mer.csv'
outh2 = open (out2,'w')
outh2.write ('#Kmer,Window,Ref,Coverage,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5\n')
for k in set (k_pool):
		Qs = ",".join (mem_window[k]['q'])
		Mis = ",".join (mem_window[k]['m'])
		Ins = ",".join (mem_window[k]['i'])
		Del = ",".join (mem_window[k]['d'])
		Cov =  ":".join (mem_window[k]['cov'])
		outh2.write (",".join ([k,Cov,Qs,Mis,Ins,Del])+'\n')
outh2.close()
os.remove (out1)
