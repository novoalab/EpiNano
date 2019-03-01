import sys
from collections import defaultdict
current = defaultdict(list)
import numpy as np
import gzip

if (len (sys.argv) !=3):
    print "USAGE:"
    print	  " python2", sys.argv[0], "per_rd_5mer_slided_current_intensity", "per_site_slided_one_kmer_oneline"
    exit()
per_rd_sliding_window_with_fast5_current_intensity = sys.argv[1]
per_site_var_slided = sys.argv[2]



#old version:
#Read,kmer,position,Quality,Mismatches,Insertions,Deletions,Reference,RefPos,Intensity
#new version:
#ReadKmer,RefKmer,Ref,RefPos,Read,ReadPos,Q,M,I,D,CurrentIntensity
if sys.argv[1].endswith ('.gz'):
		fh = gzip.open (sys.argv[1],'r')
else:
		fh = open (sys.argv[1],'r')

for line in fh:
    ary = line.strip().split(',')
    if line.startswith ('#'):
        continue
    if ary[-1] == 'NA':
        continue
    k = ary[3]+':'+ary[5]
    current[k].append (float(ary[-1]))
fh.close()


if sys.argv[2].endswith ('.gz'):
    fh = gzip.open (sys.argv[2],'r')
else:
    fh = open (sys.argv[2],'r')
header = fh.readline()
print header.strip() + ',meanCurrent,MedianCurrent,StdCurrent'

for line in fh:
    ary = line.strip().split(',')
    k = ary[2]+':'+ary[1]
    mn_intense = ''
    md_intense = ''
    std_intense = ''
    if k in current:
        intense = current[k]
        mn_intense = str (np.mean(intense))
        md_intense = str (np.median(intense))
        std_intense = str (np.std(intense))
        print ','.join ([line.strip(), mn_intense,md_intense,std_intense])
    else:
        print line.strip() + ',' + ",".join (['NA', 'NA' , 'NA'])

fh.close()
#Read,Read_Window,ReadKmer,Ref,RefKmer,Ref_Window,Q,M,I,D,CurrentIntensity
#Kmer,Window,Ref,Coverage,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5
