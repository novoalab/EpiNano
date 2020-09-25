#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import pandas as pd 
import numpy as np
from collections import defaultdict 
import sys 	

'''
Description: prepare for delta features
'''

usage = """
	python Epinano_make_delta.py <modified sample feature table> <unmodified sample feature table> <windown size of feature table>
"""

if (len (sys.argv) < 4):
	print (usage, file=sys.stderr)
	exit (0)

#~~~~~~~~~~~~~~~ arguments 
mod = sys.argv[1]
unm = sys.argv[2]
kmerLen = int(sys.argv[3])

#~~~~~~~~~~~~~~~ information to fetch 
idx = list(range(4))
#coverage 
qidx = list (range(5,5+kmerLen))
midx = list(range(5+kmerLen,5+kmerLen+kmerLen))
iidx = list (range(5+kmerLen+kmerLen,5+kmerLen+kmerLen+kmerLen))
didx = list (range(5+kmerLen+kmerLen+kmerLen,5+kmerLen+kmerLen+kmerLen+kmerLen))
#sidx = 5+kmerLen+kmerLen+kmerLen+kmerLen # samples id 

'''
only applied to 5mer 
#columns1 = list(range(4)) # Kmer,Window,Ref,Strand,
#columns2 = list(range (5,25,1))  #mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5
#columns3 = list (range(26,32,1)) # mean_I,median_I,std_I,mean_time,median_time,std_time 
26th column is the sample types information
depending on where it is put, columns number is subjected to change 
'''
modinfo = defaultdict(list)
with open (mod,'r') as fh:
	for l in fh:
		if l.startswith ('#'):
			continue 
		ary = l.rstrip().split(',')
		idx = ",".join (ary[:4])
		var = np.array (ary[qidx[0]:didx[-1]+1]).astype(float)
		#cur = np.array (ary[26:]).astype (float)
		#cur[-3:] = np.log10(cur[-3:]) # duration time transformed 
		modinfo[idx].append (var)
		#modinfo[idx].append (cur)

unminfo = defaultdict(list)
with open (unm,'r') as fh:
	for l in fh:
		if l.startswith ('#'):
			continue 
		ary = l.rstrip().split(',')
		idx = ",".join (ary[:4])
		var = np.array (ary[qidx[0]:didx[-1]+1]).astype(float)
		#cur = np.array (ary[26:]).astype (float)
		#cur[-3:] = np.log10(cur[-3:])
		unminfo[idx].append (var)
		#unminfo[idx].append (cur)

hk = '#Kmer,Window,Ref,Strand'
hq = ",".join ([f'DeltaQ{str(i+1)}' for i in range (kmerLen)])
hm = ",".join ([f'DeltaMis{str(i+1)}' for i in range (kmerLen)])
hi = ",".join ([f'DeltaIns{str(i+1)}' for i in range (kmerLen)])
hd = ",".join ([f'DeltaDel{str(i+1)}' for i in range (kmerLen)])

header = ",".join ([hk, hq, hm, hi, hd])
print (header)

delta = defaultdict(list)
for k in modinfo.keys():
	if k in unminfo: 
		var = modinfo[k][0] - unminfo[k][0]
		var = ",".join (map (str,var))
		#cur = modinfo[k][1] - unminfo[k][1]
		#cur = ",".join (map(str, cur))
		modinfo[k] = None 
		unminfo[k] = None 
		print (k,var,sep=',')
		#print (k,var,cur,sep=',')
