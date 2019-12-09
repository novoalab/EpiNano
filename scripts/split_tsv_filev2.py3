#!/usr/bin/env python3
import sys, os, math
import numpy as np
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument ("-f", "--file", required=True, help = 'input tsv file produced by sam2tsv.jar')
parser.add_argument ("-t", "--threads",type=int, help =
'number of cpus/threads that will be used to process the input file;default is 4')
parser.add_argument ('-r', '--reads', type=int, help=
'number of reads processed in batch; recommend (0,4000]; default is 1500')
parser.add_argument ('-k', '--keep_tmp_data', action='store_true')
parser.add_argument ('-v', '--verbose', action='store_true')





ary = np.array([1,2,3,4,5,6,7,8,9])
s = 2
chunks = int (math.ceil (ary.shape[0]/s))  #[1] /2
sub = np.array_split(ary,chunks)


def get_reads (tsvfile, number_of_reads_in_each_file=2000):
    reads = set()
    with open (tsvfile, 'r') as fh:
        for l in fh:
            if l.startswith ('#'):
                continue
        rd  = l.strip().splti()[0]
        reads.add(rd)
    return np.array_split (np.array (list (reads)), number_of_reads_in_each_file)





if __name__ == "__main__":
