import sys
import os
import gzip

per_read_var_slided_one_base_one_line = sys.argv[1]
mem_line = {}
#Relative_pos,Window,ReadKmer,RefKmer,Ref,RefPos,RefBase,Read,ReadPos,ReadBase,Q,M,I,D

tmp_file = sys.argv[1] + '.tmp'
tmp_out = open(tmp_file,'w')
if sys.argv[1].endswith ('.gz'):
	fh = gzip.open (sys.argv[1])
else:
	fh = open (sys.argv[1],'r')

for l in fh:
	if l.startswith('#'):
		continue
	ary =  l.strip().split (',')
	if ary[0].startswith ('0'):
		print >>tmp_out,l.strip()
tmp_out.close()
fh.close()

deletions = {}
with open (tmp_file,'r') as fh:
	for l in fh:
		if l.startswith('#'):
			continue
		ary = l.strip().split (',')
		if len (ary) < 14:
			continue
		deletion = int (ary[-1])
		if deletion > 0:
			tmp = fh.next().strip().split(',')
			new_del = int (round (deletion/2.0))
			tmp[-1] = str(new_del)
			ary[-1] = str(new_del)
			k2 = ','.join ([tmp[1],tmp[4],tmp[5],tmp[7],tmp[8]])
			k1 = ','.join ([ary[1],ary[4],ary[5],ary[7],ary[8]])
			deletions[k1] =  ",".join(ary[2:])
			deletions[k2] =  ",".join(tmp[2:])
fh.close()
os.remove (tmp_file)

if sys.argv[1].endswith ('.gz'):
	fh = gzip.open (sys.argv[1])
else:
	fh = open (sys.argv[1],'r')
for l in fh:
	if l.startswith ('#'):
		print l.strip(); continue
	ary = l.strip().split(',')
	k = ','.join ([ary[1],ary[4],ary[5],ary[7],ary[8]])
	if k in deletions:	
		print ary[0] +','+deletions[k]
	else:
		print l.strip()
fh.close()
#os.remove (sys.argv[1])
