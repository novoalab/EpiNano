import sys
from collections import defaultdict
import re

usage = ''' python ''' + sys.argv[0]+ ''' <read_length(mapped or all; trimmed or raw)>  <event_feature_tbl> \
<leading  trimming length>  \
<tail trimming length> '''

if len (sys.argv) < 4:
	print usage
	exit(0)

event = defaultdict(dict)
rd_len = {}
lead_trim = int (sys.argv[3]) # trim off first 5 bp
tail_trim = int (sys.argv[4])  # head 5 and tail 3

#inputs
fastq_len_file = sys.argv[1]
feature_event_tbl = sys.argv[2]
#sliding_win_single_rd = sys.argv[3]

with open (fastq_len_file, 'r') as fh:
    for line in fh:
        ary = line.strip().split()
        rd_len[ary[0]] = int (ary[1])

#READ,mean,stdv,m-state,move,weights,mp_state,p_state_model,model_diff,base_diff,read_pos
header = ','.join (['#Read','position','Intensity'])
print header

with open (sys.argv[2],'r') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        ary = line.strip().split(',')
        k1 = ary[0] # read id
        if k1 not in rd_len:
            continue
        positions = range (int(ary[-1]) -2 , int(ary[-1]) + 3 )
        new_positions = []
        for i in positions:
            new_positions.append (rd_len[k1] + lead_trim + tail_trim  - i - lead_trim + 1 )
        k2 = ':'.join(
            list(
            reversed (map (str,new_positions))
            ))
	print k1+','+k2+','+ary[1]
#        event[k1][k2] = ary[1]	
	
#ReadKmer,RefKmer,Ref,RefPos,RefBase,Read,ReadPos,ReadBase,Q,M,I,D  ## kmer in k lines
#ReadKmer','RefKmer','Ref','RefPos',Read','ReadPos','Q','M','I','D' ## kmer in one line

exit()

with open (sliding_win_single_rd,'r') as fh:
    lines = []
    for line in fh:
        if line.startswith('#') or line.startswith('ReadKmer'):
            continue
        ary = line.strip().split(',')
        if ary[4] not in rd_len:
            continue
        rd_km = ary[0]
        rd_pos = ':'.join (ary[5].split('|'))
	quali = ary[6]
	mis = ary[7]
	ins = ary[8]
	de = ary[9]
	rf_i = ary[2]
	rf_km = ary[1]
	rf_p = ary[3]
	#print rd_km,rd_pos,rd_q,rd_m,rd_i,rd_d
        if rd_pos in event[ary[4]]:
            print ','.join ([ary[4],rd_km,rd_pos,quali,mis,ins,de,rf_i,rf_p,rf_km,event[ary[4]][rd_pos]])
        else:
            print ','.join ([ary[4],rd_km,rd_pos,quali,mis,ins,de,rf_i,rf_p,rf_km,'NA'])
        
    
