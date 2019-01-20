import sys
from collections import defaultdict
import pandas as pd
from collections import deque

''''
different from version1
this version, aka, version2 combines deletions in reads
so as to print out complete 5mers instead of gaps in 5mers
'''

def window (seq,size=5):
    it = iter(seq)
    win = deque ((next (it,None) for _ in xrange (size)), maxlen=size)
    yield win
    append = win.append
    for e in it:
        append(e)
        yield win

# main()
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
fh = open (sys.argv[1],'r')
h = fh.readline()
row_counter = 1

out_header = ",".join(['#ReadKmer','RefKmer','Ref','RefPos',
                       'Read','ReadPos','Q','M','I','D'])
print out_header

# first line always starts with a read containing perfect match or variant other than indel;
firstLine = fh.readline()
tmp = firstLine.strip().split(',')
rd = tmp[3]
ref = tmp[0]
comb = rd+'\t'+ref
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
    if len (tmp) != 10:
	continue
    rd = tmp[3]
    ref = tmp[0]
    comb = rd+'\t'+ref
    if comb in rd_rf_mem_di:
        if int (tmp[9]) == 0 :
            rd_rf_mem_di.append(comb)
            rd_bases.append (tmp[5])
            rf_bases.append (tmp[2])
            rd_pos.append (tmp[4])
            rf_pos.append (tmp[1])
            qual.append(tmp[6])
            mis.append(tmp[7])
            ins.append(tmp[8])
            delt.append(tmp[9])
        elif int(tmp[9]) > 0:
            delt[-1] = str ( int (delt[-1]) + 1)
            rf_bases[-1] = rf_bases[-1] +'+'+tmp[2]
            #qual.append(tmp[6])
            #rd_rf_mem_di.append (comb)
            rf_pos[-1] += ';'+tmp[1]
            #rd_pos[-1] += ';'+rd_pos[-1]
            rd_bases[-1] += '-' + tmp[2]
            #qual[-1] += '-' + qual[-1]
    else:
        for i in window(range(len(rd_rf_mem_di)), 5):
            try:
                last = i[-1] + 1
            except:
                sys.stderr.write ('exception',i, rd_rf_mem_di[i], rd_pos[i],rf_pos[i])
            RD, REF =  rd_rf_mem_di[-1].split()
            rd_kmer = ''.join (rd_bases[i[0]:last])
            rf_kmer = ''.join (rf_bases[i[0]:last])
	    rf_positions = '|'.join (rf_pos[i[0]:last])
	    rd_positions = '|'.join (rd_pos[i[0]:last])
   	    qualities = '|'.join (qual[i[0]:last])
	    mismatches = '|'.join (mis[i[0]:last])
  	    insertions = '|'.join (ins[i[0]:last])
            deletions = '|'.join (delt[i[0]:last])
	    print ",".join ([rd_kmer, rf_kmer,REF,rf_positions,RD,rd_positions,qualities,mismatches,insertions,deletions])
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
    last = len (i)
    RD, REF =  rd_rf_mem_di[-1].split()
    rd_kmer = ''.join (rd_bases[i[0]:last])
    rf_kmer = ''.join (rf_bases[i[0]:last])
    for k in (i):
        print ",".join ([rd_kmer, rf_kmer, REF, rf_pos[k],rf_bases[k],
                                  RD,rd_pos[k],rd_bases[k],
                                  qual[k],mis[k],ins[k],delt[k]])

