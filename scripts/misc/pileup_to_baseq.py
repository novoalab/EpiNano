#!/usr/bin/env python
import sys
import numpy as np
from scipy.stats import mannwhitneyu

readbases = {}
qualicontainer = {}
pileup1 = sys.argv[1]
pileup2 = sys.argv[2]
with open (pileup1,'r') as p:  # native sample
    for line in p:
        line = line.rstrip("\n")
        data = line.split()
        key = ",".join(data[0:3])
        qualicontainer[key] = data[5]
p.close ()

print "#Ref,Pos,Refbase,MeanQ1,MeanQ2,MedianQ1,MedainQ2"
with open (pileup2,'r') as p: # amplified sample
    for line in p:
        line.rstrip("\n")
        data = line.split ()
        key = ",".join(data[0:3])
	if key in qualicontainer:
		ary1 = np.array ([int (ord(x) - 33 ) for x in qualicontainer[key]])
        	ary2 = np.array ([int (ord(x) - 33 ) for x in data[5]])
        	mean1,mean2 = np.mean(ary1),np.mean(ary2)
	        median1,median2 = np.median (ary1),np.median(ary2)
        	print key+','+",".join([str(x) for x in [mean1,mean2,median1,median2]])
