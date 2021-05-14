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
	python Epinano_make_delta.py <modified sample feature table> <unmodified sample feature table> <minimum coverage at sites>  <windown size of feature table>
	if not specified, slided window size is 1, aka not slided
	otherwise, windown size = kmer size
"""

if (len (sys.argv) < 5):
	print (usage, file=sys.stderr)
	exit (0)

#~~~~~~~~~~~~~~~ arguments 
mod = sys.argv[1]
unm = sys.argv[2]
cov = int (sys.argv[3])
kmerLen = int(sys.argv[4])

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
		depth = ary[4].split(':')
		middle = int(len (depth)/2)
		if int(depth[middle]) < cov:
			#print ('site', ary[:5],'is not deep enough', file=sys.stderr)
			continue
		idx = ",".join (ary[:4])
		var = np.array (ary[qidx[0]:didx[-1]+1]).astype(float)
		modinfo[idx].append (var)


hk = '#Kmer,Window,Ref,Strand'
hq = ",".join ([f'DeltaQ{str(i+1)}' for i in range (kmerLen)])
hm = ",".join ([f'DeltaMis{str(i+1)}' for i in range (kmerLen)])
hi = ",".join ([f'DeltaIns{str(i+1)}' for i in range (kmerLen)])
hd = ",".join ([f'DeltaDel{str(i+1)}' for i in range (kmerLen)])
header = ",".join ([hk, hq, hm, hi, hd])
print (header)

unminfo = defaultdict(list)
with open (unm,'r') as fh:
	for l in fh:
		if l.startswith ('#'): 
			continue 
		ary = l.rstrip().split(',')
		depth = ary[4].split(':')
		middle = int(len (depth)/2)
		if int(depth[middle]) < cov:
			#print ('site', ary[:5],'is not deep enough', file=sys.stderr)
			continue
		idx = ",".join (ary[:4])
		var = np.array (ary[qidx[0]:didx[-1]+1]).astype(float)
		if idx in modinfo:
			delta = modinfo[idx] - var 
			#delta = np.array (delta[0][:]).astype(str)
			delta = ",".join (map (str,delta[0]))
			print (idx,delta,sep=",")
		#unminfo[idx].append (var)

exit() 

'''
#delta = defaultdict(list)
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
'''
