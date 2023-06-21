#!/usr/bin/env python
import sys
import pysam
import re
import numpy as np 
from collections import defaultdict 
import multiprocessing  as mp #import Process 
from Bio import SeqIO 
import argparse

def bam_to_per_site_var(bam, refname, fafh, strands):
    '''
    bamfh: pysam.AlignemntFile bam file handle 
    refname: reference id 
    seqlen: reference sequence length dictionary 
    '''
    bamfh = pysam.AlignmentFile(bam,'rb')
    #fafh = pysam.FastaFile (fafn)
    ins, dele = {}, {} 
    for pileupcolumn in bamfh.pileup (refname, 0, bamfh.get_reference_length(refname),  flag_filter = 3844, min_mapping_quality=0, min_base_quality=0):
        pileupcolumn.set_min_base_quality (0)
        #pileupstr = pileupcolumn.get_query_sequences(add_indels=True)
        #rdnames = " ".join([r.split('-')[-1][-3:] for r in pileupcolumn.get_query_names()])
        #num_align = pileupcolumn.get_num_aligned()
        cov = num_segs = pileupcolumn.nsegments
        if cov == 0 :
            continue
        #num_reads = len(pileupcolumn.get_query_names())
        refseq = fafh.fetch (refname)
        strand = ''
        mis, mat = 0, 0 
        refname = pileupcolumn.reference_name
        refpos = pileupcolumn.reference_pos  
        refbase = refseq[refpos] if refpos is not None else 'None'
        if refpos-1 in ins: del ins[refpos-1]
        if refpos-1 in dele: del dele[refpos-1]
        ins[refpos] = ins.get(refpos,0) 
        dele[refpos] = dele.get(refpos,0) 
        base_qs = pileupcolumn.get_query_qualities()
        mean_q, median_q, std_q = '{0:.{1}f}'.format (np.mean (base_qs), 4), '{0:.{1}f}'.format(np.median(base_qs),4), '{0:.{1}f}'.format(np.std(base_qs),4)
        for pileupread in pileupcolumn.pileups:
            indel = pileupread.indel 
            qpos = pileupread.query_position
            qpos = qpos if qpos is not None else None  
            qname = pileupread.alignment.query_name.split('-')[-1][-4:]
            qseq = pileupread.alignment.query_sequence
            qbase = qseq[qpos] if qpos is not None else 'None' 
            strand = '-' if pileupread.alignment.is_reverse else '+'  
            if not strand in strands:
                continue
            if qbase != refbase and qpos:
                mis += 1 
            elif qbase == refbase and qpos:
                mat += 1
            if indel > 0:
                ins[refpos] = ins[refpos] + 1
                ins[refpos+1] = ins.get(refpos+1,0) + 1
            elif indel < 0:
                for i in range (1, abs(indel)+1):
                    dele[refpos+i] = dele.get(refpos+i,0) + 1
        mis_freq,ins_freq,dele_freq = '{0:.{1}f}'.format(mis/cov,4), '{0:.{1}f}'.format(ins[refpos]/cov,4), '{0:.{1}f}'.format(dele[refpos]/cov,4)
        print (refname, refpos + 1, refbase, strand, num_segs, mean_q, median_q, std_q, mis_freq, ins_freq, dele_freq,  sep=",")                


def main ():
    #usage = f'\n_USAGE_:\n python {sys.argv[0]} <reference file (fasta format)>   <bam file (indexed)>   <number of threads: less than or equal to the number of references>\n'
    #if len (sys.argv) != 4:
    #    print ('\nIllegal number of arguments')
    #    print (usage)
    #    sys.exit()
    
    import os  
    import glob 
    import re

    parser = argparse.ArgumentParser()
    required_args = parser.add_argument_group ('Required Arguments')
    required_args.add_argument ('-R','--reference', help='genome or transcriptome reference sequence in fasta format ') 
    required_args.add_argument ('-b', '--bam', type=str, help='bam file; and it has to be indexed')
    parser.add_argument ('-n', '--number_cpus', type=int, default=1, help='number of CPUs')
    parser.add_argument ('-T', '--type', type=str, default="t", help="reference type: either g(enome) or t(ranscriptome); default is (t)ranscriptome;")
    args=parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    fafn = args.reference #sys.argv[1]
    fafh = pysam.FastaFile (fafn)
    nreferences = fafh.nreferences
    ncpus = args.number_cpus #int(sys.argv[3])
    bam = args.bam #sys.argv[2]
    if not os.path.isfile (f'{bam}.bai'):
        print (f'\nPlease run samtools index {bam} to generate {bam}.bai \n')
        sys.exit()
    
    strands = []
    refType = args.type
    if refType == 't':
        strands = ['+']
    elif refType == 'g':
        strands = ['+', '-']
    else:
        print ('reference type has to be either t(ranscriptome) or g(enome)')
        sys.exit ()

    bamfh = pysam.AlignmentFile (bam,'rb')

    #output header
    print ('#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del')
    if ncpus  == 1 :
        ref_names = bamfh.references
        for ref in ref_names:
            bam_to_per_site_var (bam, ref, fafh, strands)
    elif ncpus > 1:
        small_files = [] 
        nreferences_per_file = int(nreferences/ncpus) 
        if nreferences_per_file == 0: # when more CPUs than number of reference sequences
            nreferences_per_file = 1
        n = 1 
        cnt = 1
        small_fn = "{}.{}".format(fafn,cnt)
        small_files.append (small_fn)
        outfh = open (small_fn, 'w')
        for rec in SeqIO.parse (fafn,'fasta'):
            if n > 1 and n % nreferences_per_file == 0: 
                outfh.close()
                cnt += 1
                small_fn = "{}.{}".format(fafn,cnt)
                small_files.append (small_fn)
                outfh = open (small_fn,'w')
            print ('>'+rec.id, rec.seq,sep="\n", file=outfh)
            n += 1
        outfh.close()
        processes = []
        for small_fn in small_files:
            pysam.faidx (small_fn)
            smallfafh = pysam.FastaFile (small_fn)
            for ref in smallfafh.references:
                processes.append (mp.Process(target=bam_to_per_site_var, args = (bam, ref, smallfafh, strands)))
        for ps in processes:
            ps.start()
        for ps in processes:
            ps.join()

    # cleaning
    if ncpus > 1:
        files = glob.glob (f'{fafn}.*')
        for f in files:
            if re.search (f'{fafn}.[0-9].fai', f):
                os.system (f'rm {f}')
            elif re.search (f'{fafn}.[0-9]', f):
                os.system (f'rm {f}')

if __name__ == '__main__':
    main ()
    #test run
    #python 2.bam2var.py -n 1 -T t -R nova_lab_epinano/test_data/make_predictions/ref.fa -b nova_lab_epinano/test_data/make_predictions/ko.bam

