#!/usr/bin/env python3
import re
import sys
from collections import defaultdict
# python oneliner.fasta
'''
usage = " USGAE:\n python motif.genome.index.py fasta_file motif_regex motif_length out \n"
if len (sys.argv) != 5:
    sys.stderr.write (usage)
    exit()
'''
import argparse
ap = argparse.ArgumentParser ()
ap.add_argument ('-m1','--forward_motif',  help = 'motif on the froward strand; e.g., RRACH would be [AG][AG]AC[ACT]')
ap.add_argument ('-m2', '--reverse_motif', help='motif comes from the reverse strand, but read on forward strand; eg., RRACH would be [AGT]GT[CT][CT]')
ap.add_argument ('-f','--fasta',required = True, help = 'fasta file')
ap.add_argument ('-l', '--motif_length',help='motif length; in the future this will be automatically computed!')
ap.add_argument ('-o','--output', required = True, help = 'output file')
args = ap.parse_args()

'''
output format:
chr, start (0-based), end (0 based), strand,  motif
'''
fasta_file = args.fasta
fmotif,rmotif = '',''
if (args.forward_motif):
    fmotif = args.forward_motif
if (args.reverse_motif):
    rmotif = args.reverse_motif
motif_len = args.motif_length
out = open (args.output,'w')

#nickNames= ['BBBBBABBBBB', 'BBBBBBBBBBB','NNANN','NNBNN','RRACH','RRBCH','BBABB','BBBBB','BSABS','BSBBS','m6a','antim6a']
#motif_len = [11,11,5,5,5,5,5,5,5,5,5,5,5]
def scan_motif (seqid, seq,  motif, motif_len):
    seq = re.sub(r'\s+','',seq.strip())
    starts=[m.start() for m in re.finditer (motif,seq)]
    inf = defaultdict(list)
    for i in starts:
        inf[seqid].append ((str(i+1),str(i+int (motif_len)), seq[i:i+int (motif_len)]  ))
    return inf

# for curlcakes on scane plus strand
seq_lst = []
seq_dic = dict()
seqid = ''
with open (args.fasta,'r') as fasta:
    for line in fasta:
        if len (line.strip()) == 0:
            continue
        if line.startswith ('>'):
            line = line.replace('>','')
            seqid = line.strip().split()[0]
            if (not seqid in seq_dic) and (len (seq_lst) > 0):
                if fmotif:
                    inf = scan_motif (seq_lst[-1], seq_dic[seq_lst[-1]], fmotif, motif_len)
                    for f in inf[seq_lst[-1]]:
                        out.write ("\t".join ([seq_lst[-1], f[0], f[1],f[2],'+\n']))
                if rmotif:
                    inf = scan_motif (seq_lst[-1], seq_dic[seq_lst[-1]], rmotif, motif_len)
                    for f in inf[seq_lst[-1]]:
                        out.write ("\t".join ([seq_lst[-1], f[0], f[1],f[2],'-\n']))
                del seq_dic[seq_lst[-1]]
            seq_lst.append (seqid)
            seq_dic[seqid] = ''
        else:
            seq = re.sub(r'\s+','',line.strip())
            seq_dic[seqid] = seq_dic[seqid] + seq
