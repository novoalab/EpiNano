#!/usr/bin/env python3
import sys
from collections import defaultdict
import pandas as pd
from collections import deque
import gzip
import os

''''
input file:
(summed up from sam2tsv results using
#REF,REF_POS,REF_BASE,READ_NAME,READ_POSITION,READ_BASE,BASE_QUALITY,MISMATCH,INSERTION,DELETION
cc6m_2244_T7_ecorv,31,A,4ea19788-8597-42c4-a7c2-3de4da5a7390,38,A,16,0,0,0
cc6m_2244_T7_ecorv,32,T,4ea19788-8597-42c4-a7c2-3de4da5a7390,39,T,3,0,0,0

this will output
one base one line format
and
one kmer one line format
'''

usage = "python " + sys.argv[0] + " per_read_var"
if len (sys.argv) <2:
    print ("USGAE: "+usage)
    exit(0)

def window (seq,size=5):
    it = iter(seq)
    win = deque ((next (it,None) for _ in range (size)), maxlen=size)
    yield win
    append = win.append
    for e in it:
        append(e)
        yield win

# main()
relative_pos = ['-2','-1','0','+1','+2']
rd_rf_mem_di = [] #defaultdict(list)
rd_bases = [] #defaultdict(str)
rf_bases = [] #defaultdict(str)
rd_pos = [] #defaultdict(str)
rf_pos = [] #defaultdict(str)
qual = [] #defaultdict(str)
mis = [] #defaultdict(str)
ins = [] #defaultdict(str)
delt = [] #defaultdict(str)

#REF", 'REF_POS','REF_BASE','READ_NAME','READ_POSITION','READ_BASE','BASE_QUALITY','MISMATCH','INSERTION','DELETION'
prefix = ''

if sys.argv[1].endswith(".gz"):
	fh = gzip.open(sys.argv[1],'r')
	prefix = '.'.join (sys.argv[1].replace('.gz','').split('.')[:-1])
else:
	fh = open (sys.argv[1],'r')
	prefix = '.'.join (sys.argv[1].split('.')[:-1])
h = fh.readline()
row_counter = 1

slided_file = prefix+'.5mer.tmp.csv'
slided_fh = open (slided_file,'w')
slided_header = '#Relative_Pos,window,ReadKmer,RefKmer,Ref,RefPos,RefBase,Read,ReadPos,ReadBase,Q,M,I,D'
slided_fh.write (slided_header+'\n')

# first line always starts with a read containing perfect match or variant other than deletion;
## slided var
firstLine = fh.readline()
tmp = firstLine.strip().split(',')
rd = tmp[3]
ref = tmp[0]
comb = rd+' '+ref
rd_rf_mem_di.append (comb)
rd_bases.append (tmp[5])
rf_bases.append (tmp[2])
rd_pos.append (tmp[4])
rf_pos.append (tmp[1])
qual.append(tmp[6])
mis.append(tmp[7])
ins.append(tmp[8])
delt.append(tmp[9])

for line in fh:
    if line.startswith('#'):
        continue
    tmp = line.strip().split(',')
    if len (tmp) < 10:
        continue
    rd = tmp[3]
    ref = tmp[0]
    comb = rd+' '+ref
    if comb in rd_rf_mem_di:
        rd_rf_mem_di.append (comb)
        rd_bases.append (tmp[5])
        rf_bases.append (tmp[2])
        rd_pos.append (tmp[4])
        rf_pos.append (tmp[1])
        qual.append(tmp[6])
        mis.append(tmp[7])
        ins.append(tmp[8])
        delt.append(tmp[9])
    else:
        for i in window(range(len(rd_rf_mem_di)), 5):
            last = None
            rd_win_pos = ''
            if i[-1] is None:  # when a read have less than window size number of positions mapped to a reference
                last = i[0] + len (i) + 1
                rd_win_pos = ':'.join(rd_pos[i[0]:last])
            else:
                last = i[-1] + 1
                rd_win_pos = ':'.join(rd_pos[i[0]:last])
		#sys.stderr.write ('exception',i, rd_rf_mem_di[i], rd_pos[i],rf_pos[i])
                #sys.stderr.write ('exception',i, rd_rf_mem_di[i], rd_pos[i],rf_pos[i])
            RD, REF =  rd_rf_mem_di[-1].split()
            rd_kmer = ''.join (rd_bases[i[0]:last])
            rf_kmer = ''.join (rf_bases[i[0]:last])
            int_ele = 0
            for k in (i):
                if isinstance (k,int):
			#sys.stderr.write (str (len (rd_rf_mem_di)) + ' rf_rd: ' + rd_rf_mem_di[0] + ' ' + rd_rf_mem_di[-1] + ' '+ ':'.join (map (str, rd_pos)) + ' ' + ':'.join (map (str, i))+'\n' )
                    slided_fh.write (",".join ([relative_pos[int_ele], rd_win_pos, rd_kmer, rf_kmer, REF, rf_pos[k],rf_bases[k],
                                 RD,rd_pos[k],rd_bases[k],
                                 qual[k],mis[k],ins[k],delt[k]]) + '\n')
                    int_ele += 1
                elif k is None: # less than 5 sites of a read mapped
                    slided_fh.write (",".join ([relative_pos[int_ele], rd_win_pos, rd_kmer, rf_kmer, REF, 'Nan,Nan',
                                RD,'Nan,Nan,Nan,Nan,Nan,Nan'])+'\n')
                    int_ele += 1
        rd_rf_mem_di = [] #defaultdict(list)
        rd_bases = [] #defaultdict(str)
        rf_bases = [] #defaultdict(str)
        rd_pos = [] #defaultdict(str)
        rf_pos = [] #defaultdict(str)
        qual = [] #defaultdict(str)
        mis = [] #defaultdict(str)
        ins = [] #defaultdict(str)
        delt = [] #defaultdict(str)
        rd_rf_mem_di.append (comb)
        rd_bases.append (tmp[5])
        rf_bases.append (tmp[2])
        rd_pos.append (tmp[4])
        rf_pos.append (tmp[1])
        qual.append(tmp[6])
        mis.append(tmp[7])
        ins.append(tmp[8])
        delt.append(tmp[9])

for i in window(range(len(rd_rf_mem_di)),5):
    last = None
    rd_win_pos = ''
    if i[-1] is None:
        last = i[0] + len(i) + 1
        rd_win_pos = ':'.join(rd_pos[i[0]:last])
    else:
        last = i[-1] + 1
        rd_win_pos = ':'.join(rd_pos[i[0]:last])
    RD, REF = rd_rf_mem_di[-1].split()
    rd_kmer = ''.join (rd_bases[i[0]:last])
    rf_kmer = ''.join (rf_bases[i[0]:last])
    int_ele = 0
    for k in (i):
        if isinstance (k, int):
            slided_fh.write (",".join ([relative_pos[int_ele],rd_win_pos,rd_kmer, rf_kmer, REF, rf_pos[k], rf_bases[k],
                                  RD,rd_pos[k],rd_bases[k],
                                  qual[k],mis[k],ins[k],delt[k]]) + '\n')
            int_ele += 1
        elif k is None:
            slided_fh.write (",".join ([relative_pos[int_ele],rd_win_pos,rd_kmer, rf_kmer, REF, 'Nan,Nan',
                                RD,'Nan,Nan,Nan,Nan,Nan,Nan']) + '\n')
            int_ele += 1

slided_fh.close()
'''
-2,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,114,G,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,80,G,6,0,0,0
-1,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,115,T,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,81,T,16,0,0,0
0,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,116,T,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,82,T,21,0,0,0
+1,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,117,T,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,83,T,23,0,0,0
+2,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,118:119:120,T+C+G,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,84,T,17,0,0,2
'''
##########################################################################
##############  sum up onebase oneline into onekmer one line format ######
##########################################################################

sum_file = prefix + '.per_rd_var.5mer.csv'
sum_out = open (sum_file,'w')

'''
input
#Relative_Pos,window,ReadKmer,RefKmer,Ref,RefPos,RefBase,Read,ReadPos,ReadBase,Q,M,I,D
-2,2:3:4:5:6,TCGAT,TCGAT,cc6m_2244_T7_ecorv,46,T,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,2,T,5,0,0,0
-1,2:3:4:5:6,TCGAT,TCGAT,cc6m_2244_T7_ecorv,47,C,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,3,C,6,0,0,0
0,2:3:4:5:6,TCGAT,TCGAT,cc6m_2244_T7_ecorv,48,G,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,4,G,12,0,0,0
+1,2:3:4:5:6,TCGAT,TCGAT,cc6m_2244_T7_ecorv,49,A,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,5,A,10,0,0,0
+2,2:3:4:5:6,TCGAT,TCGAT,cc6m_2244_T7_ecorv,50,T,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,6,T,4,0,0,0
'''
from collections import defaultdict
mem_window = defaultdict (defaultdict(list).copy)

var = slided_file
if var.endswith ('.gz'):
	f = gzip.open (var,'r')
else:
	f = open (var,'r')
header = f.readline()
new_header = ",".join (header.strip().split(',')[2:])
sum_out.write ("#Read,Read_Window,ReadKmer,Ref,RefKmer,Ref_Window,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5\n")

ary = f.readline().strip().split(',')
k =  ','.join(ary[1:5]) + ','+ary[7]
mem_window[k]['q'].append(ary[-4])
mem_window[k]['m'].append(ary[-3])
mem_window[k]['i'].append(ary[-2])
mem_window[k]['d'].append(ary[-1])
mem_window[k]['refpos'].append(ary[5])

for l in  f:
	ary = l.strip().split(',')
	k =  ','.join(ary[1:5]) + ','+ary[7]
	if k in mem_window:
		mem_window[k]['q'].append(ary[-4])
		mem_window[k]['m'].append(ary[-3])
		mem_window[k]['i'].append(ary[-2])
		mem_window[k]['d'].append(ary[-1])
		mem_window[k]['refpos'].append(ary[5])
	else:
		for kk in mem_window:
		    lst = kk.split(',')
		    h = ",".join ([lst[4],lst[0],lst[1],lst[3],lst[2]])
		    sum_out.write (','.join ([h ,':'.join (mem_window[kk]['refpos']),  ','.join (mem_window[kk]['q']),  ','.join (mem_window[kk]['m']),  ','.join (mem_window[kk]['i']),  ','.join (mem_window[kk]['d'])]) +'\n')

		del mem_window[kk]
		mem_window[k]['q'].append(ary[-4])
		mem_window[k]['m'].append(ary[-3])
		mem_window[k]['i'].append(ary[-2])
		mem_window[k]['d'].append(ary[-1])
		mem_window[k]['refpos'].append(ary[5])
for kk in mem_window.keys():
    lst = kk.split(',')
    h = ",".join ([lst[4],lst[0],lst[1],lst[3],lst[2]])
    sum_out.write  (','.join ([h ,':'.join (mem_window[kk]['refpos']),  ','.join (mem_window[kk]['q']),  ','.join (mem_window[kk]['m']),  ','.join (mem_window[kk]['i']),  ','.join (mem_window[kk]['d'])]) + '\n')
f.close()
sum_out.close()
os.remove (slided_file)
