import sys
import numpy as np
import re

header = '#READ,mean,stdv,m-state,move,weights,mp_state,p_state_model,model_diff,base_diff,read_pos'
print header

sequence = ''
mem = []

def cmp_str (str1,str2):
    '''
    compare strings of the same length
    '''
    diff = 0
    for i in range (len(str1)):
        if str1[i] == str2[i]:
            diff += 0
        else:
            diff += 1
    return diff
'''
main ()
'''

fh = open (sys.argv[1],'r')
h = fh.readline()
firstLine = fh.readline()
ary = firstLine.strip().split()
rd, mean, stdv, m_state, move, weights,p_state_model, mp_state = ary[0],ary[1],ary[2],ary[5],ary[6],ary[7],ary[8], ary[9]
m_state = re.sub ('[^agctAGCT]+','', m_state)
mp_state = re.sub ('[^agctAGCT]+','', mp_state)
max_index = np.argmax (map (float, ary[-4:]))
model_diff = cmp_str(m_state,mp_state)
bases = ['A','C','G','T']
'''
this is the order in which the probability of each base is shown in event table
so it is subjected to change
of cos, you can read from the header as well
but for now, i am lazy, so i will leave it as it is
'''
base_diff = 0 if bases[max_index] == m_state[2] else 1
mem.append ([rd, mean, stdv, m_state, move, weights, mp_state,p_state_model, model_diff, base_diff])
sequence += m_state

for line in fh:
    ary = line.strip().split()
    if ary[0] == mem[-1][0]:
        if ary[6] == '0': # skip 0 moves
            continue
        rd,mean, stdv,m_state, move,weights,p_state_model, mp_state = ary[0],ary[1],ary[2],ary[5],ary[6],ary[7],ary[8], ary[9]
        m_state = re.sub ('[^agctAGCT]+','', m_state)
        mp_state = re.sub ('[^agctAGCT]+','', mp_state)
        max_index = np.argmax (map (float, ary[-4:]))
        model_diff = cmp_str(m_state,mp_state)
        bases = ['A','C','G','T']
        base_diff = 0 if bases[max_index] == m_state[2] else 1
        mem.append ([rd,mean,stdv,m_state,move,weights,mp_state,p_state_model,model_diff,base_diff])
        sequence += m_state[-int(move):]
    else:
        firstKmer = mem.pop(0)
        kmer = sequence[0:5]
        m_state = firstKmer[2]
        print ",".join(map(str,(firstKmer))) + ',' + '3'
        pos = 3
        start = 0
        for idx in range (0,len (mem)):
            m_state = mem[idx][3]
            move = int (mem[idx][4])
            start +=  move
            kmer = sequence[start:start+5]
            print ",".join(map (str, (mem[idx]))) + ',' + str(start + pos)
            '''
            if kmer == m_state:
                print mem[idx]; print ','+str(idx+3+start_offset)
            '''
        mem = []
        ary = line.strip().split()
        if ary[6] == '0': # skip 0 moves
            continue
        rd,mean,stdv,m_state,move,weights,p_state_model, mp_state = ary[0],ary[1],ary[2],ary[5],ary[6],ary[7],ary[8], ary[9]
        m_state = re.sub ('[^agctAGCT]+','', m_state)
        mp_state = re.sub ('[^agctAGCT]+','', mp_state)
        max_index = np.argmax (map (float, ary[-4:]))
        model_diff = cmp_str(m_state,mp_state)
        bases = ['A','C','G','T']
        base_diff = 0 if bases[max_index] == m_state[2] else 1
        mem.append ([rd, mean, stdv, m_state, move, weights, mp_state,p_state_model, model_diff, base_diff])
        sequence = m_state #[-int(move):]
