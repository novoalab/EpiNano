#!/usr/bin/env python3
import sys
import gzip
from collections import defaultdict
from collections import OrderedDict
reads = {}

usage = " ". join(["python", sys.argv[0] , "<read-kmer-current-intensity> <per_read_var.5mer.csv>  <per_site_var.5mer.csv>"])

def fopen (f):
    if f.endswith ('.gz'):
        fh = gzip.open (f,'r')
    else:
        fh = open (f,'r')
    return fh

if len (sys.argv) < 3:
    print (usage)
    exit()

read_pos_intensity = sys.argv[1]
single_read_var_slided = sys.argv[2]
out = '.'.join (sys.argv[2].split('.')[:-1]) + '.current.csv'#sys.argv[3]
per_read_var_5mer_current = out
outh = open (out,'w')
##read,mean,stdv,model_state,move,read_pos
##Read,Read_Window,ReadKmer,Ref,RefKmer,Ref_Window,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5
with fopen (single_read_var_slided) as fh:
    for l in fh:
        if l.startswith ('#'):
            continue
        ary = l.strip().split(',')
        rd, win = ary[0],ary[1]
        k = rd+':'+win
        reads[k] = l.strip()

header = "#Read,Read_Window,ReadKmer,Ref,RefKmer,Ref_Window,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5,Current"
outh.write (header+'\n')
number_of_columns = len (header.split(','))
with fopen (read_pos_intensity) as fh:
    for l in fh:
        if l.startswith ('#'):
            continue
        ary = l.strip().split(',')
        rd, win = ary[0],ary[-1]
        k = rd+':'+win
        if k in reads:
            reads[k] = reads[k] + ',' + ary[1]

for k,v in reads.items():
    if len (v.split(',')) == number_of_columns:
        outh.write (v+'\n')
    else:
        outh.write (v+',NA\n')
outh.close()
##############################################################################
########### per read current intensity to per site curent intensity ##########
##############################################################################
current = defaultdict(list)
mem = OrderedDict()
import numpy as np

per_site_var_slided = sys.argv[3]
out = '.'.join (per_site_var_slided.split('.')[:-1]) + '.current.csv'
outh = open (out,'w')
# read in per site variants table and record reference location of kmers
with fopen (per_site_var_slided)as fh:
    for line in fh:
        if line.startswith ('#'):
            header = line.strip()
            continue
        ary = line.strip().split(',')
        k = ary[2]+':'+ary[1]
        mem[k] = line.strip()

with fopen (per_read_var_5mer_current) as fh:
    for line in fh:
        ary = line.strip().split(',')
        if line.startswith ('#'):
            continue
        if ary[-1] == 'NA':
            continue
        try:
            k = ary[3]+':'+ary[5]
            if k in mem:
                current[k].append(float(ary[-1]))
        except:
            sys.stderr.write (line.strip())

outh.write(header.strip() + ',meanCurrent,MedianCurrent,StdCurrent\n')

for k,v in mem.items():
    intense = current[k]
    mn_intense = str (np.mean(intense))
    md_intense = str (np.median(intense))
    std_intense = str (np.std(intense))
    if mn_intense == 'nan':
        outh.write (",".join ([v,'NA' ,"NA","NA"])+'\n')
    else:
        outh.write (','.join ([v, mn_intense,md_intense,std_intense])+'\n')
outh.close ()
