#!/usr/bin/env python3
import sys
import numpy as np
import re
import gzip

'''
positions of the middle base of the 5mers from event table in raw read is 1 based
'''

usage = "USAGE:\n python " + sys.argv[0] + ' event_tbl > event_tbl_features'
if len (sys.argv) < 2:
    print (usage)
    exit()

header = '#READ,mean,stdv,m-state,move,weights,mp_state,p_state_model,model_diff,base_diff,read_pos'
print (header)

sequence = ''
mem = []
reads = set()

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

def fopen (f):
    if f.endswith ('.gz'):
        fh=gzip.open(f,'r')
    else:
        fh = open (f,'r')
    return fh
'''
main ()
'''

fh = fopen (sys.argv[1])

'''
this is the order in which the probability of each base is shown in event table
so it is subjected to change
'''

model_diff = ''

for line in fh:
    if line.startswith ('#'):
        continue
    if re.search (r'read',line):
        continue
    ary = line.strip().split()
    rd,mean, stdv, m_state, move,weights,p_state_model, mp_state = '','','','','','','',''
    try:
        rd,mean, stdv, m_state, move,weights,p_state_model, mp_state = ary[0],ary[1],ary[2],ary[5],ary[6],ary[7],ary[8],ary[9]
    except:
        sys.stderr.write (line.strip())
    if rd.startswith('read_id'):
        continue
    if move == '0':
        continue
    if rd in reads:
        sequence += m_state[-int(move):]
        m_state = re.sub ('[^agctAGCT]+','', m_state)
        mp_state = re.sub ('[^agctAGCT]+','', mp_state)
        max_index = np.argmax (map (float, ary[-4:]))
        model_diff = cmp_str(m_state,mp_state)
        bases = ['A','C','G','T']
        base_diff = 0 if bases[max_index] == m_state[2] else 1
        mem.append ([rd,mean,stdv,m_state,move,weights,mp_state,p_state_model,model_diff,base_diff])
    else:
        if len (mem) > 0: # it is not the 1st entry of the whole file
            center = 0 # start of new read; start of new event; and start of 1sk kmer in new read
            begin = 0
            for idx in range (0,len(mem)):
                entry = mem[idx]
                m_state = entry[3]
                move = int (entry[4])
                center += move
                kmer = sequence[begin:begin+5]
                begin += move
                print (",".join(map (str,(entry))) + ',' + str(2 + center))
            mem = []
            reads.add(rd)
            rd,mean, stdv, m_state, move,weights,p_state_model,mp_state = ary[0],ary[1],ary[2],ary[5],ary[6],ary[7],ary[8],ary[9]
            seqeuence = m_state #next read first kmer
            mem.append ([rd,mean,stdv,m_state,move,weights,mp_state,p_state_model,model_diff,base_diff])
            #mem.append ([rd, mean,stdv,m_state,move,weights,p_state_model])
        else: # it is the 1st entry of the whole file with multiple reads' event information
            reads.add(rd)
            m_state = re.sub ('[^agctAGCT]+','', m_state)
            mp_state = re.sub ('[^agctAGCT]+','', mp_state)
            sequence = m_state # 1st kmer contribute wholly to sequence
            max_index = np.argmax (map (float, ary[-4:]))
            model_diff = cmp_str(m_state,mp_state)
            bases = ['A','C','G','T']
            base_diff = 0 if bases[max_index] == m_state[2] else 1
            mem.append ([rd,mean,stdv,m_state,move,weights,mp_state,p_state_model,model_diff,base_diff])

begin = 0
center = 0
for idx in range (0,len(mem)):
    entry = mem[idx]
    model_state = entry[3]
    move = int (entry[4])
    kmer = ''
    kmer = sequence[begin:begin+5]
    begin += move
    center += move
    print (",".join(map (str,(entry))) + ',' + str(2 + center))
