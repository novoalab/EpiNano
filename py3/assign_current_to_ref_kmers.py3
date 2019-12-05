#!/usr/bin/env python3
import sys
from collections import defaultdict
from collections import OrderedDict

current = defaultdict(list)
mem = OrderedDict()
import numpy as np
import gzip

'''
per read var and current format
	#Read,Read_Window,ReadKmer,Ref,RefKmer,Ref_Window,Q,M,I,D,CurrentIntensity
slided per site var format
	#Kmer,Window,Ref,Coverage,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5
'''

if (len (sys.argv) !=3):
    print (
	'''
		USAGE:
     	python3 asign_current_to_ref_kmers.py3 current_intensity  per_site_var.5mers.csv
	 '''
	 ); exit()
per_rd_current_intensity = sys.argv[1]
per_site_var_slided = sys.argv[2]

def fopen (f):
    if f.endswith('.gz'):
		fh = gzip.open (f,'r')
    else:
		fh = open (f,'r')
    return (fh)

 # read in per site variants table and record reference location of kmers
header = '#Kmer,Window,Ref,Coverage,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5'

with fopen (sys.argv[2]) as fh:
    for line in fh:
        ary = line.strip().split(',')
        k = ary[2]+':'+ary[1]
        mem[k] = line.strip()

with fopen (sys.argv[1]) as fh:
    for line in fh:
        ary = line.strip().split(',')
        if line.startswith ('#'):
            continue
        if ary[-1] == 'NA':
            continue
        try:
    	k = ary[3]+':'+ary[5]
            if k in mem:
    	    current[k].append(float(ary[-1]))
        except:
    	print line.strip()

print header.strip() + ',meanCurrent,MedianCurrent,StdCurrent'
for k,v in mem.items():
    intense = current[k]
    mn_intense = str (np.mean(intense))
    md_intense = str (np.median(intense))
    std_intense = str (np.std(intense))
    if mn_intense == 'nan':
	print ",".join ([v,'NA' ,"NA","NA"])
    else:
        print ','.join ([v, mn_intense,md_intense,std_intense])
