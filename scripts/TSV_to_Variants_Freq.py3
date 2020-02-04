#!/usr/bin/env python3
import os, sys, gzip, bz2, re
import subprocess, shutil
from multiprocessing import Pool
import argparse
from os import walk

__version__ = '2020-01-27-no-per-rd'


parser = argparse.ArgumentParser()
parser.add_argument ("-f", "--file", required=True, help = 'input tsv file produced by sam2tsv.jar')
parser.add_argument ("-t", "--threads",type=int, help =
'number of cpus/threads that will be used to process the input file;default is 4')
parser.add_argument ('-r', '--reads', type=int, help=
'number of reads processed in batch; recommend (1500,4000]; default is 1500')
parser.add_argument ('-k', '--keep_tmp_data', action='store_true')
parser.add_argument ('-p', '--per_read_stats', action='store_true', help="compute per-reads statistics")
parser.add_argument ('-v', '--verbose', action='store_true')

args = parser.parse_args()
if args.verbose:
    print(args)
    sys.exit(0)

if not args.threads:
    args.threads = 4

if not args.reads:
    args.reads = 1000
'''
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

num_of_reads_in_each_file = args.reads

tmp_dir = ''
out_prefix='small'
if os.path.isfile(args.file):
    bam_tsv = args.file
    tmp_dir = args.file + '_tmp_splitted/'
    if not os.path.exists(tmp_dir):
        os.mkdir (tmp_dir)
    else:
        #os.rmdir (tmp_dir) # only delete empty folder
        shutil.rmtree (tmp_dir)
        os.mkdir (tmp_dir)
else:
    print ('bam_tsv file '+ args.file + ' does not exist')
    exit(1)

script_dir = os.path.dirname (os.path.abspath(__file__))

smallfile = None
file_idx = 0
reads_cnt = 0
#READ_NAME      FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
reads = set()
last_reads = dict()
zero_counts = dict()

with openfile (args.file) as fh:
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
def collect_all_small_tsv_files (tsv_dir):
    tsv_files = []
    for (dirpath, dirnames,filenames) in walk(tsv_dir):
        for f in filenames:
            if re.search (r'small_\d+.tsv', f):
                tsv_file = os.path.join (dirpath, f)
                tsv_files.append(tsv_file)
    return tsv_files

small_tsv_files = collect_all_small_tsv_files (tmp_dir)

### count variants' frequencies
def tsv_to_freq (tsv): #for per site variants frequency computation
    freq_cmd = 'python3 ' + script_dir + '/per_site_var.frequency.py3 ' + tsv
    subprocess.call (freq_cmd,shell=True)
    #print (freq_cmd)

p = Pool(args.threads)
p.imap_unordered (tsv_to_freq, small_tsv_files)
p.close ()
p.join () #  Wait for all child processes to close.
#p.imap_unordered (tsv_to_freq, small_tsv_files)
# combine freq and compute variants percentage

combine_site_var_cmd = 'cat {}/small_*.freq  | python3 {}/combine_multi_per_site_var_frq.py3 > {}.per.site.var.csv'.format(tmp_dir, script_dir, args.file)
subprocess.call (combine_site_var_cmd, shell=True)
slide_win_cmd_per_site =  "python3 {}/slide_per_site_var.py3 {}.per.site.var.csv".format (script_dir, args.file)
subprocess.call (slide_win_cmd_per_site,shell=True)
## process per read variants
def tsv_to_read_var (tsv):
    cmd = 'python3 ' + script_dir +'/PER_READ.py3 ' + tsv
    subprocess.call (cmd, shell=True)

p = Pool (args.threads)
if args.per_read_stats:
    p.imap_unordered (tsv_to_read_var, small_tsv_files)
    p.close()
    p.join()
## cat per reads metrics
    cat_per_read_var_cmd = "cat {}/small_*per_rd_var.csv | awk \'NR==1 || !/^#/ \'> {}.per_rd_var.csv ".format (tmp_dir, args.file)
    subprocess.call (cat_per_read_var_cmd,shell=True)
    cat_per_read_slided_cmd = "cat {}/small_*per_rd_var.5mer.csv | awk \'NR==1 || !/^#/ \'> {}.per_rd_var.5mer.csv ".format (tmp_dir, args.file)
    subprocess.call (cat_per_read_slided_cmd,shell=True)
## delete folder
if not args.keep_tmp_data:
    shutil.rmtree (tmp_dir)
