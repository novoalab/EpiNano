import sys
import gzip
from collections import defaultdict
from collections import OrderedDict

mem_window = defaultdict (defaultdict(list).copy)
mem_window_m = defaultdict (list)
mem_window_i = defaultdict (list)
mem_window_d = defaultdict (list)
'''
input per site
#Kmer,window,Relative_pos,Ref,Ref_Pos,base,cov,q_mean,q_median,q_std,mis,ins,del
NNGTC,43:44:45:46:47,-2,Null
NNGTC,43:44:45:46:47,-1,Null
NNGTC,43:44:45:46:47,+0,cc6m_2244_T7_ecorv,45,G,1,5.000,5.000,0.000,0.000,0.000,0.000
NNGTC,43:44:45:46:47,+1,cc6m_2244_T7_ecorv,46,T,3,7.333,8.000,1.700,0.000,0.000,0.000
NNGTC,43:44:45:46:47,+2,cc6m_2244_T7_ecorv,47,C,3,5.333,6.000,1.700,0.000,0.000,0.000
NGTCG,44:45:46:47:48,-2,Null
NGTCG,44:45:46:47:48,-1,cc6m_2244_T7_ecorv,45,G,1,5.000,5.000,0.000,0.000,0.000,0.000
NGTCG,44:45:46:47:48,+0,cc6m_2244_T7_ecorv,46,T,3,7.333,8.000,1.700,0.000,0.000,0.000
NGTCG,44:45:46:47:48,+1,cc6m_2244_T7_ecorv,47,C,3,5.333,6.000,1.700,0.000,0.000,0.000
NGTCG,44:45:46:47:48,+2,cc6m_2244_T7_ecorv,48,G,3,10.000,9.000,1.414,0.000,0.000,0.000
GTCGA,45:46:47:48:49,-2,cc6m_2244_T7_ecorv,45,G,1,5.000,5.000,0.000,0.000,0.000,0.000
'''
var = sys.argv[1]

bases = defaultdict (list) 
Qs = defaultdict (list)
Mis = defaultdict (list)
Del = defaultdict (list)
Pos = defaultdict (list)
Cov = defaultdict (list)

if var.endswith ('.gz'):
	f = gzip.open (var,'r')
else:
	f = open (var,'r')

k_recorder = OrderedDict() 
for l in f:
	if l.startswith ('N'):
			continue
	elif l.startswith ('#'):
			continue
	ary = l.strip().split(',')
	if ary[0].endswith ('N'):
			continue
	ks = ary[0:4]
	del ks[2]
	k = ",".join (ks)
	k_recorder[k] = True
	mem_window[k]['q'].append(ary[7])
	mem_window[k]['m'].append(ary[10])
	mem_window[k]['i'].append(ary[11])
	mem_window[k]['d'].append(ary[12])
	mem_window[k]['cov'].append(ary[6])

f.close()
print '#Kmer,Window,Ref,Coverage,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5'
for k in k_recorder:
		Qs = ",".join (mem_window[k]['q'])
		Mis = ",".join (mem_window[k]['m'])
		Ins = ",".join (mem_window[k]['i'])
		Del = ",".join (mem_window[k]['d'])
		Cov =  ":".join (mem_window[k]['cov'])
		print ",".join ([k,Cov,Qs,Mis,Ins,Del])
