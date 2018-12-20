import sys
import numpy as np
from collections import defaultdict
from collections import OrderedDict

current = defaultdict(list)
mem = OrderedDict()
usage = '''
cat per.read.var.current.csv | python sys.argv[0] > out 
''' 

if len (sys.argv) < 2:
	print usage
	exit(1)

for line in sys.stdin:
#ReadKmer,RefKmer,Ref,RefPos,Read,ReadPos,Q,M,I,D,CurrentIntensity
	if line.startswith('#'):
		continue
	ary = line.strip().split(',')
	if ary[-1] == 'NA':
		continue
	ref,ref_pos,Intense,ref_kmer = ary[2], ary[3], ary[10],ary[1]
	k = ref+','+ref_pos+','+ref_kmer 
	current[k].append(float(Intense))
	mem[k] = True
print "#Ref,RefPos,RefKmer,MeanIntensity,MedianIntensity,Std_Intensity"
for k in mem.keys():
	print ','.join (map (str, [k,np.mean(current[k]),np.median(current[k]),np.std(current[k])]))

