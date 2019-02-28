import sys
from collections import defaultdict
from collections import OrderedDict
import numpy as np
import warnings
import re

def clear_warnings (fn):
    warnings.filterwarnings('error')
    try:
        return fn
    except RuntimeWarning:
        return 'na'

#sam2tsv file
#READ_NAME      FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
'''
reference position is 1-based
read position is 0-based

379943dd-8f04-40da-a3cf-84b6bc1a9853    0       cc6m_2595_T7_ecorv      777     A       %       873     A       M
379943dd-8f04-40da-a3cf-84b6bc1a9853    0       cc6m_2595_T7_ecorv      .       .       .       874     T       D
379943dd-8f04-40da-a3cf-84b6bc1a9853    0       cc6m_2595_T7_ecorv      .       .       .       875     C       D
379943dd-8f04-40da-a3cf-84b6bc1a9853    0       cc6m_2595_T7_ecorv      778     T       #       876     T       M

379943dd-8f04-40da-a3cf-84b6bc1a9853    0       cc6m_2595_T7_ecorv      805     A       %       902     A       M
379943dd-8f04-40da-a3cf-84b6bc1a9853    0       cc6m_2595_T7_ecorv      806     T       '       .       .       I
379943dd-8f04-40da-a3cf-84b6bc1a9853    0       cc6m_2595_T7_ecorv      810     C       ,       .       .       I
379943dd-8f04-40da-a3cf-84b6bc1a9853    0       cc6m_2595_T7_ecorv      811     G       *       903     G       M
'''
qual = defaultdict(list) # qulity scores
mis = defaultdict(int) # mismatches
ins = defaultdict(int) # insertions
dele = defaultdict(int) # deletions
cov = OrderedDict ()  # coverage
pos = defaultdict(list) # reference positions
base = {} # ref base
Q = defaultdict(list)
with open (sys.argv[1],'r') as fh:
    for line in fh:
        if line.startswith ('#'):
            continue
        ary = line.strip().split()
        if ary[-1] == 'M':
            k = (ary[2], ary[-3])
            cov[k] = cov.get(k,0) + 1
            qual[k].append (ord(ary[-4])- 33)
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
            k = cov.keys()[-1]
	    next_k = (k[0], str (int(k[1])+1))
            ins[k] = ins.get(k,0) + 0.5 #1
            ins[next_k] = ins.get(next_k,0) + 0.5 #
            #cov[k] = cov.get(k,0) + 1
            #cov[next_k] = cov.get(next_k,0) + 1
#            qual[k].append (ord(ary[-4])- 33)
#            Q[k].append(ary[-4]) # print it out and comapre it with BAQ

header = '#Ref,pos,base,cov,q_mean,q_median,q_std,mis,ins,del'
print header
for k in cov.keys():
    depth = float (cov[k])
    try:
        Mis = "%.3f" % (mis[k]/depth)
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
        print '#mean', k
    try:
        Md_q = "%.3f" % np.median(qual[k])
    except:
        print '#median', k
    try:
        sd_q = "%.3f" % np.std(qual[k])
    except:
        print '#sd', k
    #num =(mis[k],ins[k],dele[k])
    #print k,cov[k],num
    print ','.join ( map (str, [k[0],k[1],base[k],cov[k],Mn_q,Md_q,sd_q,Mis,Ins,Del])) # 'scores',",".join (map (str,qual[k])), '%'.join(Q[k])
