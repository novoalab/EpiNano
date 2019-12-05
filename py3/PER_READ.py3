#!/usr/bin/env python3
import sys
import re, os, gzip
from collections import defaultdict
import pandas as pd
from collections import deque

'''
notice that in sam2tsv output
 read position is 0-based
 reference posiiton is 1-based
input: #READ_NAME      FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
'''
if len(sys.argv) != 2:
    print ("\n usage: python PER_READ.py tsv_file \n") #ref   position\n"
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

prefix = ".".join (sys.argv[1].split('.')[:-1])
del_tmp = prefix + '.per_read.tmp.csv'
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

adjusted_file = prefix+'.per_rd_var.csv'
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

################################
### sliding per read variants
### this will output
### one base one line format
### and
### one kmer one line format
################################

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
fh = open (adjusted_file,'r')
h = fh.readline()
row_counter = 1

slided_file = prefix+'.per_read.5mer.tmp.csv'
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
