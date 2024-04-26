#!/usr/bin/env python
import sys
import re
import pysam
import numpy as np 
from collections import defaultdict 
from collections import Counter 
import multiprocessing as mp  # import Process
import argparse 
import os 
import glob
import shutil
import fileinput 
import time 


desc = """convert bam to per reference site variants frequencies.
In contrast to epinano_var1.2 and previous versions, it does not rely on sam2tsv any more.
It is also less memory demanding!
"""

epilog = '''Author: huanle.liu@crg.eu
@BCN,  10/03/2021
'''

def _rm (f):
    '''
    why can't it be defined withing main()
    https://stackoverflow.com/questions/50972934/deleting-files-with-multiprocessing-in-python
    '''
    os.remove(f)

def mkdir (folder):
    pool = mp.Pool(4)
    if os.path.exists (folder):
        if os.path.isdir (folder):
            print (folder,'exists; deleting it')
            tmp_files = glob.glob ('.csv')
            pool.map (_rm, tmp_files)
            shutil.rmtree(folder)
            os.mkdir (folder)
        else:
            print (folder, 'exists, but does not loook like a folder, please rename it to avoid being overwritten!', file=sys.stderr)
            sys.exit()
    else: 
        os.mkdir (folder)

def decimal_format (f):
    return "{0:.{1}f}".format(f,5)

def has_reads_mapped (bam) :
    return int (pysam.flagstat(bam).split('\n')[4][0]) > 0

def has_reads_mapped_to_ref (bam, refseqid):
    bamfh = pysam.AlignmentFile(bam,'rb')
    cnt = bamfh.count (contig=refseqid, start=0, end=bamfh.get_reference_length(refseqid)-1)
    return (cnt>1) #no reads mapped

def split_bam (bam):
    fwdbam = re.sub(r'bam$','fwd.bam',bam)
    revbam = re.sub(r'bam$','rev.bam',bam)
    pysam.view ("-F","16","-h","-b", "-o", fwdbam, bam, catch_stdout=False) 
    pysam.view ("-f","16","-h","-b", "-o", revbam, bam, catch_stdout=False) 
    pysam.index (fwdbam)
    pysam.index (revbam)
    return fwdbam ,revbam 

def bam_to_var (bamfn, fafn, refname, start, end, strand, outdir):
    bam = pysam.AlignmentFile (bamfn)
    fafh = pysam.FastaFile (fafn)
    if refname in fafh.references: 
        refseq = fafh.fetch (refname)
        outfn = "{}/{}".format (outdir, re.sub (r'bam$',os.path.basename(refname)+'.per.site.csv', bamfn))
        out = open (outfn, 'w')
        for pileupcolumn in bam.pileup (refname, start, end, flag_filter = 3844, min_mapping_quality=0, min_base_quality=0):
            ins = defaultdict(int) 
            pileupcolumn.set_min_base_quality (0)
            cov = pileupcolumn.get_num_aligned() if pileupcolumn.get_num_aligned() == pileupcolumn.nsegments else None
            mis, mat, dele = 0, 0 ,0
            refpos = pileupcolumn.reference_pos
            refbase = refseq[refpos].upper() #if refpos is not None else 'None'
            qq = pileupcolumn.get_query_qualities()
            mean_q, median_q, std_q = [decimal_format(x) for x in [np.mean(qq), np.median(qq), np.std(qq)]]
            print (refname, refpos, strand, refbase, cov, mean_q, median_q, std_q, end = ",", sep=",", file=out) 
            base_counts = Counter (pileupcolumn.get_query_sequences(mark_matches=False, add_indels = False))
            for base, cnt in base_counts.items():
                if len (base) == 1 and base in 'ACGTacgt':
                    if base != refbase: 
                        mis += cnt #wired '' is in 'acgt'
                    elif base == refbase:
                        mat += cnt
            base_counts =   Counter (pileupcolumn.get_query_sequences(mark_matches=True, add_indels = True))
            for base, cnt in base_counts.items ():
                if '+' in base:
                    ins[refpos] = ins.get (refpos,0) + cnt 
                if '*' == base:
                    dele = cnt
            print (decimal_format(mat/cov), decimal_format(mis/cov), decimal_format(ins[refpos]/cov), decimal_format(dele/cov),  sep=",", file=out)
            if refpos - 1 in ins: del ins[refpos-1]
    return outfn 

def multi_processes_bam2var(fafn, ref_names, bam, strand, ncpus, outdir):
    pool = mp.Pool (processes=ncpus) 
    fafh = pysam.FastaFile(fafn)
    conditions = []
    for i, j in enumerate(ref_names):
        if has_reads_mapped_to_ref(bam, j) and j in fafh.references:
            conditions.append ([bam, fafn, j,0,len(fafh.fetch(j)), strand, outdir])
    files = pool.starmap (bam_to_var,conditions)
    finaloutfn = re.sub(r'bam$','',os.path.basename(bam)) +'per.site.csv'
    finalout = open (finaloutfn,'w')
    #finalout2 = open (finaloutfn2, 'w')
    print ("#Ref,pos,strand,base,cov,q_mean,q_median,q_std,mat,mis,ins,del",file=finalout)
    try:
        for l in fileinput.input(files):
            print (l.rstrip(), file=finalout)
    except:
        numberlines = len (list (fileinput.input(files)))
        print (l, bam,'has', numberlines, 'lines', file=sys.stderr)

    
def main ():
    usage = "%(prog)s -v"
    parser = argparse.ArgumentParser(description=desc, epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-b', '--bam', required=True, type=str, help='bam file with *bai index')
    parser.add_argument('-r', '--reference', required=True, type=str, help='reference file with samtools faidx index')
    parser.add_argument('-c', '--cpus', type=int, default=4, help='number of CPUs to use [%(default)s]')
    args = parser.parse_args()

    fafn = args.reference
    bam = args.bam 
    ncpus = args.cpus
    start = time.time()
    bamfh = pysam.AlignmentFile (bam, 'rb')
    ref_names = bamfh.references
    fwdbam, revbam = split_bam(bam)
    outdir = bam.replace('.bam','') + '.tmp'
    mkdir (outdir)
    strands = ['+', '-']
    for i, bam in enumerate ([fwdbam, revbam]):
        if has_reads_mapped (bam):
            multi_processes_bam2var(fafn, ref_names, bam, strands[i], ncpus,  outdir)
    end = time.time()
    print ('analysis took {} seconds'.format(end - start ))
    files = glob.glob("{}/*.csv".format (outdir))
    pool = mp.Pool (ncpus)
    pool.map (_rm, files)
    shutil.rmtree (outdir)

    pool.map (_rm, [fwdbam, revbam, fwdbam+'.bai', revbam+'.bai'])
if __name__ == "__main__":
    main () 
