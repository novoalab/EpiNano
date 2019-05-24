import sys
from collections import defaultdict
import re
import gzip

usage = ''' python ''' + sys.argv[0]+ ''' <read_length(mapped or all; trimmed or raw)>  <event_feature_tbl> \
<leading  trimming length>  \
<tail trimming length> '''

if len (sys.argv) < 5:
	print usage
	exit(0)

def fopen (f) :
    if f.endswith ('.gz'):
	fh = gzip.open (f,'r')
    else:
	fh = open (f,'r')
    return fh

event = defaultdict(dict)
rd_len = {}
lead_trim = int (sys.argv[3]) # trim off first 5 bp
tail_trim = int (sys.argv[4])  # head 5 and tail 3

#inputs
fastq_len_file = sys.argv[1]
feature_event_tbl = sys.argv[2]
fh1 = fopen (fastq_len_file)
if (True):
    for line in fh1:
        ary = line.strip().split()
        rd_len[ary[0]] = int (ary[1])
fh1.close()
header = ','.join (['#Read','position','Intensity','Kmer_3_to_5'])
print header

fh = fopen (sys.argv[2])
#READ,mean,stdv,m-state,move,weights,mp_state,p_state_model,model_diff,base_diff,read_pos
if (True):
    for line in fh:
        if line.startswith('#'):
            continue
        ary = line.strip().split(',')  
        k1 = ary[0] # read id
        if k1 not in rd_len:
            continue
        positions = range (int(ary[-1]) -2 , int(ary[-1]) + 3 )  # position of middle base in 5mers is 1based
        new_positions = []
        for i in positions:
            new_positions.append (rd_len[k1] + lead_trim + tail_trim  - i - lead_trim + 1 )
        k2 = ':'.join(
            list(
            reversed (map (str,new_positions))
            ))
	print k1+','+k2+','+ary[1] +','+ary[3]
fh.close()
exit()
