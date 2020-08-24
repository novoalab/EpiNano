#!/usr/bin/env python

import sys, gzip

plus_strand = dict()

fh = None
if sys.argv[1].endswith('.gz'):
	fh = gzip.open (sys.argv[1], 'rt')
else:
	fh = open (sys.argv[1])

out1 = gzip.open (sys.argv[2]+'.forward_strand.gz','wt')
out2 = gzip.open (sys.argv[2]+'.reverse_strand.gz','wt')

for l in fh:
	ary = l.rstrip().split()
	plus_strand[ary[0]] = True	
fh.close()


fh = None
if sys.argv[2].endswith('.gz'):
        fh = gzip.open (sys.argv[2], 'rt')
else:
        fh = open (sys.argv[2])

for l in fh:
	if l.startswith('contig'):
		print (l.rstrip(),file=out1)
		print (l.rstrip(),file=out2)
		continue 
	ary = l.rstrip().split()
	if ary[3] in plus_strand:
		print (l.rstrip(),file=out1)
	else:
		print (l.rstrip(),file=out2)
fh.close()
	
	
	


