import sys
import pandas as pd
import numpy as np
import gzip

from itertools import groupby, islice, tee
from operator import itemgetter

def openfile (fn):
    if fn.endswith ('.gz'):
        return gzip.open (fn,'rt')
    return open (fn,'r')

def consecutive_subseq(iterable, length):
    for _, consec_run in groupby(enumerate(iterable), lambda x: x[0] - x[1]):
        k_wise = tee(map(itemgetter(1), consec_run), length)
        for n, it in enumerate(k_wise):
            next(islice(it, n, n), None) # consume n items from it
        yield from zip(*k_wise)

def fetch_consecutive_sites (lines):
    pos = sorted ([int(x.split(',')[1]) for x in lines])
    contig = list(set ([x.split(',')[0] for x in  lines]))
    if len (contig) == 1 and len (*consecutive_subseq (pos, 5) ) == 5:
        win = f'{pos[0]}-{pos[-1]}'
        mer = "".join ([x.split(',')[2] for x in lines])
        quals = ",".join ([x.split(',')[5] for x in lines])
        mis = ",".join ([x.split(',')[8] for x in lines])
        ins = ",".join ([x.split(',')[9] for x in lines])
        dels = ",".join ([x.split(',')[10] for x in lines])
        covs = ":".join ([x.split(',')[4] for x in lines])
        strand = lines[0].split(',')[3]
        return mer,win,list(contig)[0],strand,covs,quals,mis,ins,dels
    else:
        return None

''''
input file:
#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del
output file:

'''
def main():
    usage = f'''  _USAGE_:\n\n    python {sys.argv[0]} <input> <window size>
    e.g.: python Slide_Variants.py ko.per.site.csv 5\n 
    input is the output of Epinano_Variants.py;
    input has header:#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del
    input should be sorted on contig and position using: sort -t',' -k1,1 -k2,2n

    output will be directed to a file named {input}.per.site.5mer.csv

    Please post relevant issue on epinano github repository.
    Author: Huanle Liu 
    Email: elzedliu@gmail.com
    '''


    if len (sys.argv) != 3:
        print ('\n Illegal number of arguments!')
        print (usage, file = sys.stderr)
        sys.exit()
    
    inp = sys.argv[1]
    winsize = int(sys.argv[2])
    
    if 'per.site.csv' in inp:
        i = inp.index ('per.site.csv')
        out = open (f'{inp[:i]}per.site.5mer.csv', 'w')
    elif 'csv' in inp:
        i = inp.index ('csv')
        out = open (f'{inp[:i]}per.site.5mer.csv', 'w')
    else:   
        out = open (f'{sys.argv[1]}.per.site.5mer.csv', 'w')
    err = open (f'{sys.argv[1]}.non-consecutive-sites', 'w')
    print ('#Kmer,Window,Ref,Strand,Coverage,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5', file=out)

    with openfile (sys.argv[1]) as fh:
        fh.readline() # skip header
        row1, row2, row3, row4, row5 = [fh.readline().strip() for _ in range (winsize)] #fh.readline().strip(), fh.readline().strip(),\
#                fh.readline().strip(),fh.readline().strip(),fh.readline().strip()
        lines = [row1, row2, row3, row4, row5]
        inf = fetch_consecutive_sites (lines)
        if inf:
            print (",".join (inf), file=out)
        else:
            print ("\n".join (lines), file=err)

        for l in fh:
            del lines[0]
            lines.append (l.strip())
            try:
                inf = fetch_consecutive_sites (lines)
                if inf:
                    print (",".join(inf), file=out)
            except:
                print ("\n".join (lines),file=err)
                continue
    
    #del lines[0]
    #if len (lines) == 5:
    #    inf = fetch_consecutive_sites (lines)
    #    if inf:
    #        print (",".join (inf), file=out)
    out.close()
    err.close()
if __name__ == '__main__':
    main ()
