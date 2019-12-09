#!/usr/bin/env python3
import sys
import re
from collections import defaultdict
import os
import gzip
'''
notice that in sam2tsv output
 read position is 0-based
 reference posiiton is 1-based
input: #READ_NAME      FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
'''
if len(sys.argv) != 2:
    print ("\n usage: python per_read_var.py tsv_file \n") #ref   position\n"
    exit(0)

sam2tsv = sys.argv[1]
rdnames = []
qualities = {}
dels = defaultdict(int)
mis = {}
ins = {}
match = {}
ref_pos = {}

k = ''
next_k = ''

if sam2tsv.endswith (".gz"):
    fh = gzip.open (sam2tsv)
else:
    fh = open (sam2tsv)
#with open (sam2tsv,'r') as fh:
for line in fh:
    if re.match ('\s+',line):
        continue
    if re.match('#',line):
        continue
    if re.match (':',line):
        continue
    ary = line.strip().split()
    if not re.match (r'[MID]',ary[-1]):
        continue
    if len (ary) != 9:
        continue
    if ary[6].startswith('-'):
        continue
    if re.match (r'[HS]',ary[-1]):
        continue
#READ_NAME      FLAG    CHROM   READ_POS  BASE  QUAL   REF_POS REF   OP
    if ary[6] != '.':
        ary[6] = str(int(ary[6]))  # ref pos is 1-based
    ref = ary[2]
    ref_pos = ary[6]
    ref_base = ary[7]
    rd = ary[0]
    rd_pos = '' #ary[3]
    if ary[-1] != 'D':
        rd_pos = str (int (ary[3]) + 1)#turn read_pos into 1-based
        rd_base = ary[4]
        k = ','.join ([ref,ref_pos,rd,rd_pos]) #include reference  to account for multi-mappings
        qualities[k] = ord(ary[5]) - 33
        dels[k] = dels.get(k,0) + 0
        rdnames.append(k)
        if ary[-1].upper() == 'M' and ary[4] == ary[7]: # or re.match ('I', ary[-1].upper()):
            mis[k] = '0'
            ins[k] = '0'
        elif ary[-1].upper() == 'M' and ary[4] != ary[7]: # or re.match ('I', ary[-1].upper()):
            mis[k] = '1'
            ins[k] = '0'
        elif ary[-1].upper() == 'I':
            ins[k] = '1'
            mis[k] = '0'
    elif ary[-1].upper() == 'D':
        dels[k]= dels.get(k,0) + 1
        mis[k] = '0'
        ins[k] = '0'
    else:
        continue


prefix = '.'.join (sys.argv[1].split('.')[:-1])
del_tmp = prefix + '.per_read_var.tmp.csv'
tmp_fh = open (del_tmp,'w')

if sam2tsv.endswith (".gz"):
    fh = gzip.open (sam2tsv)
else:
    fh = open (sam2tsv)

for l in fh:
       if l.startswith ('#'):
           continue
       if re.match ('\s+',l):
           continue
       if re.match('#',l):
           continue
       if re.match (':',l):
           continue
       ary = l.strip().split()
       if len (ary) != 9:
           continue
       if ary[6].startswith('-') :
           continue
       ref = ary[2]
       ref_pos = ary[6]
       ref_base = ary[7]
       rd = ary[0]
       rd_pos = '' #ary[3]
       if not re.match (r'[DMI]',ary[-1]):
           continue
       if all ([ary[-1] != 'D' , ary[-1] != 'S' , ary[-1] != 'H']):
           rd_pos = str (int (ary[3]) + 1)
           rd_base = ary[4]
           k = ','.join ([ref,ref_pos,rd,rd_pos])
           inf = k.split(',')
           tmp_fh.write (','.join ([inf[0],inf[1],ref_base,inf[2],inf[3],rd_base,str (qualities[k]),mis[k],ins[k],str(dels[k])])+'\n')
tmp_fh.close()

adjusted_file = prefix+'_per_rd_var.csv'
adj_fh = open (adjusted_file,'w')
adj_fh.write ( ",".join (["#REF",'REF_POS','REF_BASE','READ_NAME','READ_POSITION','READ_BASE','BASE_QUALITY','MISMATCH','INSERTION','DELETION']) + '\n')
last_del = defaultdict(float)
line_num = 0

with open (del_tmp,'r') as FH:
    for l in FH:
        ary = l.strip().split(',')
        half = float(ary[-1]) / 2
        ary[-1] = last_del.get (line_num,0)  + half
        line_num += 1
        last_del[line_num] = half
        adj_fh.write ( ','.join(map (str, ary)) + '\n')
adj_fh.close()
os.remove (del_tmp)
