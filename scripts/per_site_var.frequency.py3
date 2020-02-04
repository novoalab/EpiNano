#!/usr/bin/env python3
import sys, re, gzip
from collections import defaultdict
from collections import OrderedDict
import numpy as np
import warnings


__version__ = '0.1-2020-01-14'
__Author__ = 'Huanle.liu@crg.eu'

def clear_warnings (fn):
    warnings.filterwarnings('error')
    try:
        return fn
    except RuntimeWarning:
        return None

def openfile(f):
    if f.endswith ('.gz'):
        fh = gzip.open (f,'rt')
    elif f.endswith ('bz') or f.endswith ('bz2'):
        fh = bz2.open(f,'rt')
    else:
        fh = open(f,'rt')
    return fh

#sam2tsv file
#READ_NAME     FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
'''
reference position is 1-based
read position is 0-based
'''

usage = "python " + sys.argv[0] + " tsv_file_converted_from_bam/sam_file \n"
if (len (sys.argv) != 2):
    print ("\n __USAGE__:\n "+usage); exit()

prefix= '.'.join (sys.argv[1].split('.')[:-1])
outh = open (prefix+'.freq','w')

qual = defaultdict(list) # qulity scores
mis = defaultdict(int) # mismatches
mat = defaultdict (int) #matches
ins = defaultdict(int) # insertions
dele = defaultdict(int) # deletions
cov = OrderedDict ()  # coverage
#ins_q = defaultdict(lambda: defaultdict (list)) # store quality scores of inserted read bases
ins_q = defaultdict(list)
aln_mem = []  #read, ref, refpos; only store last entry not matching insertion
pos = defaultdict(list) # reference positions
base = {} # ref base
Q = defaultdict(list)
with openfile (sys.argv[1]) as fh:
    for line in fh:
        if line.startswith ('#'):
            continue
        ary = line.strip().split()
        if ary[-1] == 'M':
            k = (ary[2], int (ary[-3])) #
            cov[k] = cov.get(k,0) + 1
            aln_mem = []
            aln_mem.append((ary[0],ary[2],int(ary[-3])))
            qual[k].append (ord(ary[-4])- 33)
            Q[k].append(ary[-4])
            base[k] = ary[-2]
            if (ary[-2] != ary[4]):
                mis[k] += 1
            else:
                mat[k] += 1
        if ary[-1] == 'D':
            k = (ary[2], int(ary[-3]))
            cov[k] = cov.get(k,0) + 1
            aln_mem = []
            aln_mem.append((ary[0],ary[2],int(ary[-3])))
            base[k] = ary[-2]
            dele[k] = dele.get(k,0) + 1
        if ary[-1] == 'I':
            last_k = aln_mem[-1][1],aln_mem[-1][2] # last alignment with match/mismatch/del
            next_k = (ary[2], last_k[1] + 1)
            if last_k[0] != ary[2]:
                sys.stderr.write (line.strip())
            ins_k_up = (ary[0], ary[2], last_k[1])
            ins_k_down = (ary[0], ary[2], last_k[1] + 1)
            if not (ins_k_down) in ins_q:
                ins[next_k] = ins.get(next_k,0) + 1
                #cov[k] = cov.get(k,0) + 0
                #sys.stderr.write (str(ins_k_up) + "has insertion\n")
                #sys.stderr.write (str(ins_k_down) + "has insertion\n")
                ins_q[ins_k_down].append(ord(ary[-4]) - 33)
            if not (ins_k_up) in ins_q:
                ins[last_k] = ins.get(last_k,0) + 1
                #cov[k] = cov.get(k,0) + 0
                #ins_q[ins_k_up].append(ord(ary[-4]) - 33)
                ins_q[ins_k_up].append(ord(ary[-4]) - 33)

#header = '#Ref,pos,base,cov,mis,ins,del,q_sum'
for k in cov.keys():
    depth = float (cov[k])
    Mis = mis[k]
    Mat = mat[k]
    Del = dele[k]
    num_ins = 0
    if k in ins:
        num_ins = ins[k] # len (ins_q[k].keys())
    q_lst = [0]
    if k in qual:
        q_lst = qual[k]
    inf = map (str, [k[0], k[1], base[k], depth,  Mat, Mis, num_ins, Del, ':'.join (map (str, q_lst))])
    outh.write (",".join ( inf) + '\n')

    #print ','.join ( map (str, [k[0],k[1],base[k],cov[k],Mn_q,Md_q,sd_q,Mis,Ins,Del])) # 'scores',",".join (map (str,qual[k])), '%'.join(Q[k])
