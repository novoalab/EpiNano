#!/usr/bin/env python3
import os, sys, gzip, bz2
import subprocess, shutil

'''
Usage:
    python3 <split_bam_tsv_to_small_files.py> big.tsv number-of-reads-in-eac-small-file
    by defautlt 1000 reads will be allocated into each small file
    let huanle (elzedliu@gmail.com) know if you run into any type of errors
'''

def openfile(f):
    if f.endswith ('.gz'):
        fh = gzip.open (f,'rt')
    elif f.endswith ('bz') or f.endswith ('bz2'):
        fh = bz2.open(f,'rt')
    else:
        fh = open(f,'rt')
    return fh

num_of_reads_in_each_file = 1000
if len (sys.argv) < 2:
    print (usage)
    exit(1)

tmp_dir = ''
out_prefix='small'
if os.path.isfile(sys.argv[1]):
    bam_tsv = sys.argv[1]
    tmp_dir = sys.argv[1] + '_tmp_splitted/'
    if not os.path.exists(tmp_dir):
        os.mkdir (tmp_dir)
    else:
        #os.rmdir (tmp_dir) # only delete empty folder
        shutil.rmtree (tmp_dir)
        os.mkdir (tmp_dir)
else:
    print ('bam_tsv file '+ sys.argv[1] + ' does not exist')
    exit(1)

if len (sys.argv) == 3:
    try:
        num_of_reads_in_each_file = int (sys.argv[2])
    except:
        print ('please input an integer number of reads in each small file')
        raise

scritp_dir = os.path.dirname (os.path.abspath(__file__))

cat_cmd = 'cat ' + tmp_dir + 'small_*.tsv' + '|python3 ' + scritp_dir+'/combine_multi_per_site_var_frq.py3 > ' + sys.argv[1] + 'per.site.var.freq.csv'
print (cat_cmd)
exit(0)

smallfile = None
file_idx = 0
reads_cnt = 0
#READ_NAME      FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
reads = set()
last_reads = dict()
zero_counts = dict()

with openfile (sys.argv[1]) as fh:
    for l in fh:
        if l.startswith('#'):
            continue
        rd = l.split()[0]
        reads.add(rd)
        reads_cnt = len(reads)
        if reads_cnt % num_of_reads_in_each_file == 1:
            if rd not in last_reads:
                file_idx += 1
                last_reads[rd] = True
            small_filename = tmp_dir+'/'+out_prefix+'_{}.tsv'.format(file_idx)
            smallfile = open (small_filename,'a')
        smallfile.write (l)
smallfile.close()

#####  count variant frequencies
#subprocess.call (cat_cmd, shell=True)
