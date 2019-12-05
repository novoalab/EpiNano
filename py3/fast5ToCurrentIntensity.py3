#!/usr/bin/env python3
import sys, warnings,os
import numpy as np
import re
import gzip
import h5py, ont_fast5_api
from ont_fast5_api.fast5_interface import get_fast5_file

warnings.filterwarnings("ignore",category=DeprecationWarning)

usage = '''
USAGE:
    python fast5ToCurrentIntensity.py fast5file/fast5folder

fast5file is generated with albacore 2.1.7+, which should contain event table.
fast5file can be either single-read or multi-reads fast5.
E.g, event table in fast5 files produced by albacore 2.1.7 contain the following columns:
mean,stdv,start,length,model_state,move,weights,p_model_state,mp_state,p_A,p_C,p_G,p_T
Ideally, the event table should at least contain the following columns:
model_state, mean, stdv and move.
'''

if len (sys.argv) < 2:
    print ("\n",usage,'\n')
    exit(0)


def events_tbl_from_f5 (hdf5):
    '''
    for single fast5 file, it is hdf5 file open with h5py
    for multi- fast5 file, it is hdf5['read_id'] for single read
    '''
    if 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq' not in hdf5:
        return None, None
    elif 'Analyses/Basecall_1D_000/BaseCalled_template/Events' in hdf5:
        header = hdf5['Analyses/Basecall_1D_000/BaseCalled_template/Events'].value.dtype.names
        header='\t'.join (list (header))
        events = hdf5['Analyses/Basecall_1D_000/BaseCalled_template/Events'].value
        return (header, events)

def write_event_tbl (rd, header, events, fh):
    fh.write('#Read\t'+header+'\n')
    for i in events:
            i = str(i).replace("(","")
            i = str(i).replace(",","\t")
            pat1,pat2,pat3 = """b'""", """'""", """)"""
            i = i.replace (pat1,'')
            i = i.replace (pat2,'')
            i = i.replace (pat3,'')
            tbl_fh.write (str (rd) +'\t'+i+'\n')
            tbl_fh.write (str (rd) +'\t'+i+'\n')

#### main
f5files = []
dir_input = None
if os.path.isdir(sys.argv[1]):
    dir_input = True
    sys.stderr.write("{} is a folder, will extract current intensity information from all fast5 files within it\n".format(sys.argv[1]) )
    for root, dirs, files in os.walk (sys.argv[1]):
        for file in files:
            if file.endswith ('.fast5'):
                f5files.append (file)
elif os.path.isfile(sys.argv[1]):
    f5files.append (sys.argv[1])

basename = os.path.dirname (sys.argv[1])
tbl = os.path.abspath(sys.argv[1])  + '.tbl'
tbl_fh = open (tbl,'w')

for f5file in f5files:
    if dir_input:
        f5file = os.path.join (sys.argv[1], f5file)
    hdf5 = h5py.File(f5file,'r')
    f5 = get_fast5_file (f5file,mode='r')
    reads_in_f5 = f5.get_read_ids()
    number_of_reads =  len( reads_in_f5)
    # fast5 groups: <KeysViewHDF5 ['Analyses', 'PreviousReadInfo', 'Raw', 'UniqueGlobalKey']>
    if number_of_reads == 1:
        header, events = events_tbl_from_f5 (hdf5)
        read_id = reads_in_f5[0]
        if header is not None  and events is not None :
            write_event_tbl (read_id, header, events, tbl_fh)
        else:
            sys.stderr.write("no event table for read {} \n" .format(read_id))
    elif number_of_reads > 1:
        for rd in hdf5.keys():
            header, events = events_tbl_from_f5 (hdf5[rd])
            if header is not None and events is not None:
                rd = rd.replace('read_','')
                write_event_tbl (rd, header, events, tbl_fh)
            else:
                sys.stderr.write("no event table for read {} \n" .format(rd))
tbl_fh.close()

## extract current intensity values from event tables
current_intensity_file = os.path.abspath(sys.argv[1]) + '.current.csv'
current_intensity_fh = open (current_intensity_file, 'w')
#header = '#READ,mean,stdv,m-state,move,weights,mp_state,p_state_model,model_diff,base_diff,read_pos'
header = '#read,mean,stdv,model_state,move,Read_Window\n'
current_intensity_fh.write (header)

sequence = ''
mem = []
reads = set()

def cmp_str (str1,str2):
    '''
    compare strings of the same length
    between model-state and mp_state
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

fh = fopen (tbl)
model_diff = ''
idx1,idx2,idx3,idx4 = '','','',''
ary=fh.readline().strip().split()
idx1, idx2, idx3, idx4  = ary.index('mean'), ary.index('move'), ary.index('stdv'), ary.index ('model_state')
for line in fh:
    if line.startswith ('#'): #header lines of event table
        continue
    ary = line.strip().split()
    #rd,mean,stdv,m_state,move,weights,p_state_model,mp_state = '','','','','','','',''
    rd, mean,move,std,model_state = ary[0], '','','',''
    mean,move,stdv,model_state = ary[idx1], ary[idx2],ary[idx3], ary[idx4]
    if move == '0':
        continue
    if rd in reads:
        sequence += model_state[-int(move):]
        model_state = re.sub ('[^agctAGCT]+','', model_state)
        mem.append ([rd,mean,stdv,model_state,move])
    else:
        if len (mem) > 0: # it is not the 1st entry of the whole file
            center = 0 # start of new read; start of new event; and start of 1sk kmer in new read
            begin = 0
            for i in range (0,len(mem)):
                entry = mem[i]
                model_state = entry[3]
                move = int (entry[4])
                center += move
                kmer = sequence[begin:begin+5]
                begin += move
                win = ':'.join (map (str, list(range(begin,begin+5))))
                current_intensity_fh.write(",".join(map (str,(entry))) + ',' + win + '\n')
            mem = []
            reads.add(rd)
            seqeuence = model_state #next read first kmer
            mem.append ([rd,mean,stdv,model_state,move])
        else: # it is the 1st entry of the whole file with multiple reads' event information
            reads.add(rd)
            model_state = re.sub ('[^agctAGCT]+','', model_state)
            mem.append ([rd,mean,stdv,model_state,move])
            sequence = model_state # 1st kmer contribute wholly to sequence


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
    win = ':'.join (map (str, list(range(begin,begin+5))))
    #current_intensity_fh.write (",".join(map (str,(entry))) + ',' + str(2 + center)+','+win+'\n') 2019 10 16
    current_intensity_fh.write (",".join(map (str,(entry))) + ','+win+'\n')
#os.remove (tbl)
current_intensity_fh.close()
os.remove (tbl)
