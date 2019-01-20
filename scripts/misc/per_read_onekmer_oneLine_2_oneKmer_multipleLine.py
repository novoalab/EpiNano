import sys
import gzip
import re
import numpy as np

print "#Read_kmer,Relative_pos,Read,ReadPos,ReadWin,Q,M,I,D,CurrentIntensity"

relative_pos = map (str,range (-2,3))

if sys.argv[1].endswith('.gz'):
	fh =  gzip.open (sys.argv[1],'r')
else:
	fh = open (sys.argv[1],'r')

for l in fh:
	ary = l.strip().split(',')
	rd_kmer,rd,rd_win,Q,M,I,D,Current = (ary[0],ary[4],ary[5],ary[6],ary[7],ary[8],ary[9],ary[10])
	if re.search (r'[^AGCT]',rd_kmer):
		continue
	else:
	#	print ",".join ([rd_kmer, rd,rd_win,Q,M,I,D,Current])
		rd_pos = rd_win.split ('|')
		rd_q = Q.split('|')
		rd_M = M.split('|')
		rd_I  = I.split('|')
		rd_D = D.split('|')
		rd_C = Current
		for i in range (5):
			print ",".join ([rd_kmer, relative_pos[i],rd, rd_pos[i],rd_win, rd_q[i], rd_M[i], rd_I[i], rd_D[i], rd_C])
#			print ",".join ([rd,rd_pos,rd_q,rd_M,rd_I,Current])
fh.close()
		
#ReadKmer,RefKmer,Ref,RefPos,Read,ReadPos,Q,M,I,D,CurrentIntensity






