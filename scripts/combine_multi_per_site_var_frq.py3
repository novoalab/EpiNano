#!/usr/bin/env python3
import  fileinput as fin
from collections import defaultdict
from collections import OrderedDict
import numpy as np
import sys 

#cc6m_2244_T7_ecorv,40,T,6.0,0,0,0,11:7:3:11:11:7
#ref,ref_pos,ref_base,depth,mis,ins,del,Qs
mem = defaultdict(lambda: defaultdict(list))
Var=defaultdict (list)
Qual = defaultdict (list) 

ks = OrderedDict()
for l in fin.input():
    ary = l.strip().split(',')
    k = ','.join(ary[:3])
    ks[k] = True
    c,M,m,i,d = map (float,ary[3:8])
    c_m_i_d = np.array ([c,m,i,d])
    q_lst = [ float (x) for x in ary[8].split(':')]
    Var[k] = Var.get(k,np.array([0])) + c_m_i_d 
    Qual[k] = Qual.get(k,[]) + q_lst 
   
print ( '#Ref,pos,base,cov,q_mean,q_median,q_std,mis,ins,del')
for k in ks:
    cov = Var[k][0]
    q_lst = Qual[k]
    var_freq = Var[k][1:]/cov 
    var_freq = ",".join (var_freq.astype (str) )
    print ("{},{},{},{},{},{}".format (k,cov,'%0.5f'%np.mean(q_lst),'%0.5f'%np.median(q_lst),'%0.5f'%np.std(q_lst), var_freq)) 

