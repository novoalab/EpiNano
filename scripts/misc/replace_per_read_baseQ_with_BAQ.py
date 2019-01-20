#!/usr/bin/env python
import sys
import numpy as np
from collections import defaultdict
import gzip

qualicontainer = defaultdict  (list)

pileup = sys.argv[1]
if pileup.endswith ('.gz'):
	p = gzip.open(pileup,'r')
else:
	p = open (pileup,'r')

mem_rd = {} # keep record of reads from pileup
for line in p:
	line = line.rstrip("\n")
	data = line.split()
	key = ":".join(data[0:2])
	#baq = ",".join (map (str, [int(ord(x) -33) for x in data[5]]))
	baq = map (str, [int(ord(x) -33) for x in data[5]])
	qualicontainer[key] = baq

'''
#kmer,position,relative_pos,Ref,Ref_pos,Ref_base,Depth,mean_q,median_q,std_q,mis,ins,del,mean_current,median_current,std_current
'''
single_read_var = sys.argv[2]
if single_read_var.endswith (".gz"):
	fh = gzip.open  (sys.argv[2],'r')
else:
	fh = open (sys.argv[2],'r')

for l in fh:
	if l.startswith('#'):
		print l.strip()
		continue
	ary = l.strip().split(',')
	ref,ref_pos = (ary[0],ary[1])
	base_q = ary[-4]
	k = ref+':'+ref_pos
	if k in qualicontainer:
		baq = qualicontainer[k]
		try:
			ary[-4] = baq[0] 
			baq.pop(0)
			qualicontainer[k]  = baq
		except:
			sys.stderr.write('seems do not have enough baq ==> more reads in tsv than from pileup\n')
	print ",".join(ary)
		
