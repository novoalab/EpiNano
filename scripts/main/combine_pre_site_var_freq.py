import  fileinput as fin
from collections import defaultdict
from collections import OrderedDict
import numpy as np

#cc6m_2244_T7_ecorv,40,T,6.0,0,0,0,11:7:3:11:11:7
#ref,ref_pos,ref_base,depth,mis,ins,del,Qs
mem = defaultdict(lambda: defaultdict(list))
ks = OrderedDict() 
for l in fin.input():
    ary = l.strip().split(',')
    k = ','.join(ary[:3])
    ks[k] = True
    c,m,i,d = (ary[3:7])
    c_m_i_d = np.array (map(float,[c,m,i,d]))
    q_lst = map (float, ary[7].split(':'))
    mem[k]['var'] = mem[k].get ('var',np.array((0,))) + c_m_i_d
    mem[k]['q'] = mem[k].get ('q',[]) + q_lst

print '#Ref,pos,base,cov,q_mean,q_median,q_std,mis,ins,del'
for k in ks:
    cov = mem[k]['var'][0]
    q_lst = mem[k]['q']
    print ",".join ([k, str(cov),  ",".join (map (str, ['%0.3f' % np.mean(q_lst),'%0.3f'%np.median(q_lst),'%0.3f'%np.std(q_lst)])), ",".join (map (str, ['%0.3f'% x for x in list (mem[k]['var'])[1:]/cov] ))] )
    
