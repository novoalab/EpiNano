#!/usr/bin/env python3
import sys, re, os, gzip
from collections import deque
from collections import defaultdict
from collections import OrderedDict
import numpy as np
import warnings
import re

#sam2tsv file
#READ_NAME      FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
'''
reference position is 1-based
read position is 0-based
'''
prefix = ".".join (sys.argv[1].split('.')[:-1])
var_file = prefix + '.per_site.var.csv'
outh = open (var_file,'w')

qual = defaultdict(list) # qulity scores
mis = defaultdict(int) # mismatches
ins = defaultdict(int) # insertions
dele = defaultdict(int) # deletions
cov = OrderedDict ()  # coverage
ins_q = defaultdict(lambda: defaultdict (list)) # store quality scores of inserted read bases
pos = defaultdict(list) # reference positions
base = {} # ref base
Q = defaultdict(list)
with open (sys.argv[1],'r') as fh:
    for line in fh:
        if line.startswith ('#'):
            continue
        ary = line.strip().split()
        if ary[-1] == 'M':
            k = (ary[2], ary[-3]) #
            cov[k] = cov.get(k,0) + 1
            qual[k].append (ord(ary[-4]) - 33)
            Q[k].append(ary[-4])
            base[k] = ary[-2]
            if (ary[-2] != ary[4]):
                mis[k] +=1
        if ary[-1] == 'D':
            k = (ary[2], ary[-3])
            cov[k] = cov.get(k,0) + 1
            base[k] = ary[-2]
            dele[k] = dele.get(k,0) + 1
        if ary[-1] == 'I':
            k = list (cov.keys())[-1]
            next_k = (ary[2], str(int(k[-1]) + 1))
            if ary[0] not in ins_q[k]:
                ins[k] = ins.get(k,0) + 0.5
                ins[next_k] = ins.get(next_k,0) + 0.5
            ins_q[k][ary[0]].append(ord(ary[-4]) - 33)
            #ins_q[next_k][ary[0]].append(ord(ary[-4]) - 33)
        #Q[k].append(ary[-4])

        #qual[k].append (ord(ary[-4])- 33)

header = '#Ref,pos,base,cov,q_mean,q_median,q_std,mis,ins,del'
outh.write( header+'\n')

for k in cov.keys():
    depth = float (cov[k])
    Mis = mis[k]
    Del = dele[k]
    q_lst = []
#    inf = map (str, [k[0], k[1], base[k], depth, Mis, num_ins, Del, np.sum(q_lst)])
    #out.write (",".join ( inf) + '\n')
    try:
        Mis = "%.3f" % (Mis/depth)
    except:
        Mis = 0
    try:
        Ins = "%.3f" % (ins[k]/depth)
    except:
        Ins = 0
    try:
        Del = "%.3f" % (dele[k]/depth)
    except:
        Del = 0
    try:
        Mn_q = "%.3f" % np.mean(qual[k])
    except:
        sys.stderr.write ('#mean' + k + '\n')
    try:
        Md_q = "%.3f" % np.median(qual[k])
    except:
        sys.stderr.write ('#median', k)
    try:
        sd_q = "%.3f" % np.std(qual[k])
    except:
        sys.stderr.write ('#sd', k)
    #num =(mis[k],ins[k],dele[k])
    #print k,cov[k],num
    outh.write (','.join ( map (str, [k[0],k[1],base[k],cov[k],Mn_q,Md_q,sd_q,Mis,Ins,Del]))+'\n')
outh.close()

###############################
### slide per site variants ###
###############################
mem = {}
contents = OrderedDict()

win_size = 5
dist = int (win_size) // 2 + 1
fh = open (var_file,'r')
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

out1 = prefix +'.per_site.slided.tmp'
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
out2 = prefix + '.per_site.slided.csv'
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
