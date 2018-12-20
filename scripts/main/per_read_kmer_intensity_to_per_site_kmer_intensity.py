import sys
from collections import defaultdict
current = defaultdict(list)
import numpy as np

per_rd_sliding_window_with_fast5_current_intensity = sys.argv[1]
#Read,kmer,position,Quality,Mismatches,Insertions,Deletions,Reference,RefPos,Intensity

with open (sys.argv[1], 'r') as fh:
    for line in fh:
        ary = line.strip().split(',')
        if line.startswith ('#'):
            continue
        if ary[-1] == 'NA':
            continue
        k = ary[-3]+':'+ary[-2]
        current[k].append (float(ary[-1]))

per_site_var_slided = sys.argv[2]

header = '#Kmer,position,relative_pos,Ref,Ref_pos,Ref_base,Depth,mean_q,median_q,std_q,mis,ins,del,mean_current,median_current,std_current'
print header
with open (sys.argv[2],'r') as fh:
    for line in fh:
        ary = line.strip().split(',')
        k = ary[3]+':'+ary[1]
        mn_intense = ''
        md_intense = ''
        std_intense = ''
        if k in current:
            intense = current[k]
            mn_intense = str (np.mean(intense))
            md_intense = str (np.median(intense))
            std_intense = str (np.std(intense))
        print ','.join ([line.strip(), mn_intense,md_intense,std_intense])

