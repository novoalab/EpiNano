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

this will output one base one line format 
and
one kmer one line format


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
	prefix = sys.argv[1].replace('.gz','')
else:
	fh = open (sys.argv[1],'r')
	prefix = sys.argv[1]
h = fh.readline()
row_counter = 1

slided_file = prefix+'.slided.onebase.oneline'
slided_fh = open (slided_file,'w')
slided_header = '#Relative_Pos,window,ReadKmer,RefKmer,Ref,RefPos,RefBase,Read,ReadPos,ReadBase,Q,M,I,D'
print >>slided_fh, slided_header

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
	if int (tmp[9]) == 0 :
	    rd_rf_mem_di.append (comb) 
            rd_bases.append (tmp[5])
            rf_bases.append (tmp[2])
            rd_pos.append (tmp[4])
            rf_pos.append (tmp[1])
            qual.append(tmp[6])
            mis.append(tmp[7])
            ins.append(tmp[8])
            delt.append(tmp[9])
        elif int(tmp[9]) > 0:
            delt[-1] = str(int(delt[-1]) + 1)
            rf_bases[-1] += '+' + tmp[2]
            rf_pos[-1] += ':' + tmp[1]    	
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
                	print >>slided_fh,  ",".join ([relative_pos[int_ele], rd_win_pos, rd_kmer, rf_kmer, REF, rf_pos[k],rf_bases[k],
                                 RD,rd_pos[k],rd_bases[k],
                                 qual[k],mis[k],ins[k],delt[k]])
			int_ele += 1
		elif k is None: # less than 5 sites of a read mapped
			print >>slided_fh,  ",".join ([relative_pos[int_ele], rd_win_pos, rd_kmer, rf_kmer, REF, 'Nan,Nan',
                                RD,'Nan,Nan,Nan,Nan,Nan,Nan'])
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
#    sys.stderr.write(str(len(rd_rf_mem_di)) + ' ')
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
	        print >> slided_fh, ",".join ([relative_pos[int_ele],rd_win_pos,rd_kmer, rf_kmer, REF, rf_pos[k], rf_bases[k],
                                  RD,rd_pos[k],rd_bases[k],
                                  qual[k],mis[k],ins[k],delt[k]])
		int_ele += 1
	elif k is None:
		print >> slided_fh, ",".join ([relative_pos[int_ele],rd_win_pos,rd_kmer, rf_kmer, REF, 'Nan,Nan',
                                RD,'Nan,Nan,Nan,Nan,Nan,Nan'])
		int_ele += 1

slided_fh.close()
### adjust deletion counts
#### choose only centre lines
'''
-2,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,114,G,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,80,G,6,0,0,0
-1,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,115,T,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,81,T,16,0,0,0
0,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,116,T,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,82,T,21,0,0,0
+1,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,117,T,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,83,T,23,0,0,0
+2,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,118:119:120,T+C+G,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,84,T,17,0,0,2
'''
tmp_file = slided_file + '.tmp'
tmp_out = open(tmp_file,'w')

with open (slided_file,'r') as fh:
	for l in fh:
		if l.startswith('#'):
			continue
		ary =  l.strip().split (',')
		if ary[0].startswith ('0'):
			print >>tmp_out,l.strip()
tmp_out.close()
## allocate del to adjacent position
deletions = {}
err_file = tmp_file +'.err'
err_fh = open (err_file,'w')
with open (tmp_file,'r') as fh:
	for l in fh:
		if l.startswith('#'):
			continue
		ary = l.strip().split (',')
		if len (ary) < 14:
			continue
		deletion = '' 
		try:
			deletion = int (ary[-1])
		except:
			print >> err_fh,l.strip() 
			continue
		if deletion > 0:
			tmp = fh.next().strip().split(',')
			new_del = int (round (deletion/2.0))
			tmp[-1] = str(new_del)
			ary[-1] = str(new_del)
			k2 = ','.join (tmp[4:9])
			k1 = ','.join (ary[4:9])
			deletions[k1] =  ",".join(ary[4:])
			deletions[k2] =  ",".join(tmp[4:])
err_fh.close()
fh.close()
os.remove (tmp_file)

adjusted_slide_file = slided_file + '.del.adjusted'
out_fh = open (adjusted_slide_file,'w')
with open (slided_file,'r') as fh:
	for l in fh:
		if l.startswith ('#'):
			print >>out_fh, l.strip(); continue
		ary = l.strip().split(',')
		k = ','.join (ary[4:9])
		if k in deletions:	
#			sys.stderr.write(k+'\n')
			print >>out_fh, ",".join (ary[0:4]) +','+ deletions[k]
		else:
			print >>out_fh, l.strip()
out_fh.close()
os.remove (slided_file)

##########################################################################
##############  sum up onebase oneline into onekmer one line format ######
##########################################################################

sum_file = adjusted_slide_file + '.summed.oneKmer_oneLine'
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

var = adjusted_slide_file
if var.endswith ('.gz'):
	f = gzip.open (var,'r')
else:
	f = open (var,'r')
header = f.readline()
new_header = ",".join (header.strip().split(',')[2:])
print >>sum_out, "#Read_Window,ReadKmer,RefKmer,Ref,Read,Qualities,Mismatches,Insertions,Deletions"

ary = f.readline().strip().split(',')
k =  ','.join(ary[1:5]) + ','+ary[7]
mem_window[k]['q'].append(ary[-4])
mem_window[k]['m'].append(ary[-3])
mem_window[k]['i'].append(ary[-2])
mem_window[k]['d'].append(ary[-1])

for l in  f:
	ary = l.strip().split(',')
	k =  ','.join(ary[1:5]) + ','+ary[7]
	if k in mem_window:
		mem_window[k]['q'].append(ary[-4])
		mem_window[k]['m'].append(ary[-3])
		mem_window[k]['i'].append(ary[-2])
		mem_window[k]['d'].append(ary[-1])
	else:
		for kk in mem_window:
			print >>sum_out, ','.join ([kk ,'|'.join (mem_window[kk]['q']),  '|'.join (mem_window[kk]['m']),  '|'.join (mem_window[kk]['i']),  '|'.join (mem_window[kk]['d'])])
		del mem_window[kk]
		mem_window[k]['q'].append(ary[-4])
		mem_window[k]['m'].append(ary[-3])
		mem_window[k]['i'].append(ary[-2])
		mem_window[k]['d'].append(ary[-1])
for kk in mem_window.keys():
	print >>sum_out, ','.join ([kk ,'|'.join (mem_window[kk]['q']),  '|'.join (mem_window[kk]['m']),  '|'.join (mem_window[kk]['i']),  '|'.join (mem_window[kk]['d'])])		
f.close()
sum_out.close()

