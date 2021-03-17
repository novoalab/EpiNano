#!/usr/bin/env python
from collections import defaultdict
from collections import OrderedDict
import numpy as np
import pysam, sys, re, os, re, gzip, bz2, datetime 
import argparse as ap
from collections import  deque
from itertools import repeat 
import shutil,datetime 
import fileinput 
import subprocess 

__version__ = '0.1-2020-04-04'
__Author__ = 'Huanle.liu@crg.eu'

def openfile(f):
    if f.endswith ('.gz'):
        fh = gzip.open (f,'rt')
    elif f.endswith ('bz') or f.endswith ('bz2'):
        fh = bz2.open(f,'rt')
    else:
        fh = open(f,'rt')
    return fh
	
def now ():
	return datetime.datetime.now().strftime('%D:%H:%M:%S')
	
def window (seq,size=5):
    it = iter(seq)
    win = deque ((next (it,None) for _ in range (size)), maxlen=size)
    yield win
    append = win.append
    for e in it:
        append(e)
        yield win
	
def reference_from_bam (bam):
	# bam is pysam.AlignmentFile read bam 
	bamfh = pysam.AlignmentFile (bam, 'rb')
	return set (bamfh.header.references)
def filt_bam_with_pysam (bam_in):
	# filt out non-primary alignments 
	# bam_in is a pysam.AlignmentFile bam file hand 
	# bam_out output bam fileanme 
	bamfh = pysam.AlignmentFile(bam_in,'rb')
	bam_out = re.sub (r'.bam$','',bam_in) + '.filt.bam'
	outfh = pysam.AlignmentFile(bam_out,'wb',header=bamfh.header)
	bai_out = bam_out + '.bai'
	for read in bamfh.fetch():
		cond1 = read.mapping_quality < 1
		cond2 = read.is_duplicate
		cond3 = read.is_qcfail
		cond4 = read.is_secondary
		cond5 = read.is_supplementary
		cond6 = read.is_unmapped
		if any ([cond1,cond2,cond3,cond4,cond5,cond6]):
			continue
		outfh.write (read)
	outfh.close()
	pysam.index (bam_out,bai_out)
	return bam_out, bai_out
	
def filt_bam (bam):
	'''
	with csamtools 
	filt out bad alignemnts 
	'''
	out = bam.replace('bam','filt.bam') if bam.endswith ('bam') else bam +'.filt.bam'
	pysam.view ('-F', '3844', '-h','-b','-o',out, bam, catch_stdout=False)
	pysam.index (out,out+'.bai')
	return out, out+'.bai'
	
def split_bam (bam,refid):
	'''
	split bam file on individual reference sequence 
	'''
	newbam = re.sub(r'.bam$','',bam)+'.{}.bam'.format(refid)
	pysam.view (bam, refid, '-h','-b', '-o', newbam, catch_stdout = False)
	pysam.index (newbam,newbam+'.bai')
	return newbam, newbam+'.bai'

def clean_soft_hard_clippings (ref_query_pair):
    '''
    soft-clipped: bases in 5' and 3' of the read are NOT part of the alignment.
    hard-clipped: bases in 5' and 3' of the read are NOT part of the alignment AND those bases have been removed from the read sequence in the BAM file. The 'real' sequence length would be length(SEQ)+ count-of-hard-clipped-bases
    :param ref_query_pair: list of tuples, each tuple contains read_pos, ref_pos, ref_base;
    '''
    for x,y in enumerate (ref_query_pair):
        if y[1] == None and y[2] == None:
            continue
        else:
            return ref_query_pair[x:]
            break

def variant_typing (ref_query_pair_tuple):
    '''
    :param ref_query_pair_tuple:
    '''
    if isinstance (ref_query_pair_tuple[0], int) and ref_query_pair_tuple[1] is None and ref_query_pair_tuple[2] == None:
        return 'I'
    elif ref_query_pair_tuple[0]  is None and (ref_query_pair_tuple[1], int)  and ref_query_pair_tuple[2] in 'AGCT':
        return 'D'
    elif isinstance(ref_query_pair_tuple[0], int)  and isinstance( ref_query_pair_tuple[1], int)  and ref_query_pair_tuple[2] in 'agct':
        return 'M'
    elif isinstance (ref_query_pair_tuple[0], int)  and isinstance(ref_query_pair_tuple[1], int)  and ref_query_pair_tuple[2] in 'AGCT':
        return 'm'  #

def bam_to_tsv (bam):
	bamfh = pysam.AlignmentFile(bam,'rb')
	out_tsv_fh = open (bam + '.tsv', 'w')
	header = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("#READ_NAME","FLAG","CHROM","READ_POS","BASE", "QUAL","REF_POS","REF","OP",'STRAND')
	out_tsv_fh.write(header)
	for read in bamfh.fetch():
		o1, o2, o3  = read.query_name, read.flag, read.reference_name
		query_seq = read.query_sequence
		pairs = read.get_aligned_pairs(with_seq=True)
		pairs = clean_soft_hard_clippings(pairs)
		pairs = clean_soft_hard_clippings(pairs[::-1])
		pairs = pairs[::-1]
		strand = '-' if read.is_reverse else '+'
		op =''
		for p in pairs:
			try:
				o9 = variant_typing(p)
				op = o9
			except:
				sys.stderr.write ("{}\t{}\t{} is problematic\n".format (read.reference_name, read.query_name, p) )
				exit()
			if op in ['D']:
				o4, o5, o6 = '.', '.', '.'
				o7, o8 = p[1] + 1, p[2]
			elif op in ['I'] :
				o4,o5,o6 = p[0],query_seq[int(p[0])],read.query_qualities[p[0]]
				o7,o8 = '.','.'
			else:
				o4,o5,o6, o7, o8= p[0], query_seq[int(p[0])].upper(), read.query_qualities [p[0]], int (p[1]) + 1, p[2].upper()
			out_tsv_fh.write ("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(o1,o2,o3,o4,o5,o6,o7,o8,op,strand))
	out_tsv_fh.close()
	return (bam+'.tsv')

def spot_empty_tsv (tsv):
	ary = [] 
	cnt = 0 
	with open (tsv,'r')  as fh:
		for l in fh:
			if cnt <2:
				ary.append (l)
			else: 
				break 
			cnt += 1
	return True if len (ary)>1 else False 		

def split_tsv (tsv, tmp_dir, number_of_reads_in_each_file=3000):
	output_prefix = 'small'
	small_files = set ()
	smallfile = None
	file_idx = 0
	last_seen = ''
	new_start= ''
	reads_cnt = 0
	zero_counts = dict()
	small_filename = "{}/{}_{}.tsv".format (tmp_dir,output_prefix,file_idx)
	small_files.add (small_filename)               
	smallfile = open (small_filename,'a')
	with openfile (tsv) as fh:
		for l in fh:
			if l.startswith ('#'):
				continue
			rd = l.split()[0]
			if rd != last_seen:
				last_seen = rd
				reads_cnt += 1
				#print (rd, reads_cnt)
			if all ([reads_cnt > number_of_reads_in_each_file, reads_cnt % number_of_reads_in_each_file == 1, rd != new_start]):
				smallfile.close()
				new_start = rd
				file_idx += 1
				small_filename = "{}/{}_{}.tsv".format (tmp_dir,output_prefix,file_idx)
				smallfile = open (small_filename,'a')
				small_files.add (small_filename)
			smallfile.write(l)
	smallfile.close()
	sys.stderr.write ("{} reads splitted to {} files\n".format (reads_cnt, len(small_files)))
	return (small_files)

def split_tsv_for_per_site_var_freq(tsv, q, number_threads,  num_reads_per_chunk=4000):
	head = next(tsv)
	firstline = next (tsv) 
	current_rd = firstline.split()[0]

	rd_cnt = 1
	idx = 0
	chunk_out = [] # open ("CHUNK_{}.txt".format(idx),'w')
	chunk_out.append(firstline)
	try:
		for line in tsv:
				rd = line.split()[0]
				if current_rd != rd:
						rd_cnt += 1
						current_rd = rd
						if ((rd_cnt-1) % num_reads_per_chunk == 0 and rd_cnt >= num_reads_per_chunk):
								q.put ((idx, chunk_out)) #.close()
								idx += 1
								chunk_out = [] #open ("CHUNK_{}.txt".format(idx),'w')
				chunk_out.append(line)
		q.put ((idx, chunk_out))
	except:
		raise 
		sys.stderr.write("split tsv file on reads failed\n")
	finally:
		for _ in range(number_threads ):
			q.put(None)

def split_tsv_for_per_site_var_freq_1(tsv, q, number_threads,  number_of_reads_in_each_file=2000):
	'''only computing per read features need small tsv files to be kept
	input is a generator
	fh = iter ([]) 
	if isinstance (tsv,str):  #intput is tsv text file 
		fh = openfile (tsv)  
	elif isinstance (tsv,subprocess.Popen): #return from subproces.Popen 
		fh = tsv.stdout
	elif hasattr(tsv, 'read'):
		fh = tsv 
	'''
	small_chunk = []
	idx = 0
	last_seen = ''
	reads_cnt = 0
	new_start = ''
	try:
		if True:
			for l in tsv:
				if l and l.startswith ('#'):
					continue
				rd = l.split()[0]
				if rd != last_seen:
					last_seen = rd
					reads_cnt += 1
				if  (reads_cnt > number_of_reads_in_each_file and reads_cnt % number_of_reads_in_each_file == 1): #, rd != new_start]):
					idx += 1
					new_start = rd
					q.put ((idx,small_chunk))
					small_chunk=[]
					small_chunk.append(l.strip())
				small_chunk.append (l.strip())
		q.put((idx,small_chunk))
	except:
		raise 
		sys.stderr.write("split tsv file on reads failed\n")
	finally:
		for _ in range(number_threads ):
			q.put(None)
		
def tsv_to_freq_multiprocessing_without_manager (tsv_reads_chunk_q, out_dir):
	'''
	produced with sam2tsv.jar with strand information added 
	read read-flags	reference	read-pos	read-base	read-qual	ref-pos	ref-base		cigar-op		strand
	a3194184-d809-42dc-9fa1-dfb497d2ed6a    0       cc6m_2244_T7_ecorv      0       C       #       438     G       S       +
	''' 
	for idx, tsv_small_chunk in iter (tsv_reads_chunk_q.get,None):
		#sys.stderr.write("idx-{}\tq-{}\n".format(idx, tsv_small_chunk))
		filename = "{}/small_{}.freq".format(out_dir, idx)
		#small_freq_file_list.append(filename)
		outh = open (filename,'w')
		mis = defaultdict(int) # mismatches
		mat = defaultdict (int) #matches
		ins = defaultdict(int) # insertions
		dele = defaultdict(int) # deletions
		cov = OrderedDict ()  # coverage
		ins_q = defaultdict(list)
		aln_mem = []  #read, ref, refpos; only store last entry not matching insertion
		pos = defaultdict(list) # reference positions
		base = {} # ref base
		qual = defaultdict(list)
		#READ_NAME     FLAG    CHROM   READ_POS  BASE   QUAL  REF_POS REF  OP 	STRAND 
		#read read-flags	reference	read-pos	read-base	read-qual	ref-pos	ref-base		cigar-op		strand
		for line in tsv_small_chunk:
			if line.startswith ('#'):
				continue
			ary = line.strip().split()
			if ary[-2] in ['M','m']:
				k = (ary[2], int (ary[-4]), ary[-1]) #
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				qual[k].append (ord(ary[-5])-33)
				base[k] = ary[-3].upper()
				if (ary[-3] != ary[4]):
					mis[k] += 1
				else:
					mat[k] += 1
			if ary[-2] == 'D':
				k = (ary[2], int(ary[-4]), ary[-1])
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				base[k] = ary[-3].upper() 
				dele[k] = dele.get(k,0) + 1
			if ary[-2] == 'I':
				last_k = aln_mem[-1][1],aln_mem[-1][2],aln_mem[-1][3] # last alignment with match/mismatch/del
				next_k = (ary[2], last_k[1] + 1,last_k[2])
				if last_k[0] != ary[2]:
					sys.stderr.write (line.strip())
				ins_k_up = (ary[0], ary[2], last_k[1],last_k[2])
				ins_k_down = (ary[0], ary[2], last_k[1] + 1,last_k[2])
				if (ins_k_down) not in ins_q:
					ins[next_k] = ins.get(next_k,0) + 1
					ins_q[ins_k_down].append(ord(ary[-5])-33)
				if (ins_k_up) not in ins_q:
					ins[last_k] = ins.get(last_k,0) + 1
					ins_q[ins_k_up].append(ord(ary[-5])-33)
		header = '#Ref,pos,base,cov,mat,mis,ins,del,qual,strand\n'
		#cc6m_2244_T7_ecorv,7,A,1.0,1,0,0,0,15,+
		outh.write(header)
		for k in cov.keys():
			depth = cov.get (k,0) 
			Mis = mis.get (k,0) 
			Mat = mat.get (k,0)
			Del = dele.get (k,0)
			q_lst = qual.get (k,[0])
			try:
				q_lst = ':'.join (map (str, q_lst)) +':' #qual[k]  #becos of dask sum function 
				num_ins = ins.get (k,0)
				inf = "{},{},{},{},{},{},{},{},{},{}\n".format (k[0], k[1], base[k], depth, Mat, Mis, num_ins, Del, q_lst, k[2])
				outh.write (inf)
			except:
				sys.stderr.write ("file {} {} does not work\n".format (tsv,k))
				
				
def tsv_to_freq_multiprocessing_with_manager (tsv_reads_chunk_q, out_dir):
	'''
	mutliprocessing 
	produced with sam2tsv.jar with strand information added 
	read read-flags	reference	read-pos	read-base	read-qual	ref-pos	ref-base		cigar-op		strand
	a3194184-d809-42dc-9fa1-dfb497d2ed6a    0       cc6m_2244_T7_ecorv      0       C       #       438     G       S       +
	''' 
	for idx, tsv_small_chunk in iter (tsv_reads_chunk_q.get, None):
		filename = "{}/small_{}.freq".format(out_dir, idx)
		#tsv_file = open ("{}/small_{}.tsv".format(out_dir, idx),'w')
		outh = open (filename,'w')
		mis = defaultdict(int) # mismatches
		mat = defaultdict (int) #matches
		ins = defaultdict(int) # insertions
		dele = defaultdict(int) # deletions
		cov = OrderedDict ()  # coverage
		ins_q = defaultdict(list)
		aln_mem = []  #read, ref, refpos; only store last entry not matching insertion
		pos = defaultdict(list) # reference positions
		base = {} # ref base
		qual = defaultdict(list)
		#READ_NAME     FLAG    CHROM   READ_POS  BASE   QUAL  REF_POS REF  OP 	STRAND 
		#read read-flags	reference	read-pos	read-base	read-qual	ref-pos	ref-base		cigar-op		strand
		for line in tsv_small_chunk:
			if line.startswith ('#'):
				continue
			#tsv_file.write (line+'\n')
			ary = line.strip().split()
			#sys.stdout.write(line+'\n')
			if ary[-2] in ['M','m']:
				k = (ary[2], int (ary[-4]), ary[-1]) #
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				qual[k].append (ord(ary[-5])-33)
				base[k] = ary[-3].upper()
				if (ary[-3] != ary[4]):
					mis[k] += 1
				else:
					mat[k] += 1
			if ary[-2] == 'D':
				k = (ary[2], int(ary[-4]), ary[-1])
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				base[k] = ary[-3].upper() 
				dele[k] = dele.get(k,0) + 1
			if ary[-2] == 'I':
				#print (aln_mem)
				last_k = aln_mem[-1][1],aln_mem[-1][2],aln_mem[-1][3] # last alignment with match/mismatch/del
				#last_k = list (cov.keys())[-1]
				next_k = (ary[2], last_k[1] + 1,last_k[2])
				if last_k[0] != ary[2]:
					pass
					#sys.stderr.write (line.strip())
				ins_k_up = (ary[0], ary[2], last_k[1],last_k[2])
				ins_k_down = (ary[0], ary[2], last_k[1] + 1,last_k[2])
				if (ins_k_down) not in ins_q:
					ins[next_k] = ins.get(next_k,0) + 1
					ins_q[ins_k_down].append(ord(ary[-5])-33)
				if (ins_k_up) not in ins_q:
					ins[last_k] = ins.get(last_k,0) + 1
					ins_q[ins_k_up].append(ord(ary[-5])-33)
		header = '#Ref,pos,base,cov,mat,mis,ins,del,qual,strand\n'
		#cc6m_2244_T7_ecorv,7,A,1.0,1,0,0,0,15,+
		outh.write(header)
		for k in cov.keys():
			depth = cov.get (k,0) 
			Mis = mis.get (k,0) 
			Mat = mat.get (k,0)
			Del = dele.get (k,0)
			q_lst = qual.get (k,[0])
			try:
				q_lst = ':'.join (map (str, q_lst))+':'  # dataframe sum  
				num_ins = ins.get (k,0)
				inf = "{},{},{},{},{},{},{},{},{},{}\n".format (k[0], k[1], base[k], depth, Mat, Mis, num_ins, Del, q_lst, k[2])
				outh.write (inf)
			except:
				sys.stderr.write ("file {} {} does not work\n".format (tsv,k))			

def tsv_to_freq (tsv):
	'''
	single thread
	produced with sam2tsv.jar && with strand information added 
	read read-flags	reference	read-pos	read-base	read-qual	ref-pos	ref-base		cigar-op		strand
	a3194184-d809-42dc-9fa1-dfb497d2ed6a    0       cc6m_2244_T7_ecorv      0       C       #       438     G       S       +
	'''
	out = re.sub(r'.tsv$','',tsv) + '.freq' 
	outh = open (out,'w')
	mis = defaultdict(int) # mismatches
	mat = defaultdict (int) # matches
	ins = defaultdict(int) # insertions
	dele = defaultdict(int) # deletions
	cov = OrderedDict ()  # coverage
	ins_q = defaultdict(list)
	aln_mem = []  #read, ref, refpos; only store last entry not matching insertion
	pos = defaultdict(list) # reference positions
	base = {} # ref base
	qual = defaultdict(list)
#READ_NAME     FLAG    CHROM   READ_POS  BASE   QUAL  REF_POS REF  OP 	STRAND 
#read read-flags	reference	read-pos	read-base	read-qual	ref-pos	ref-base		cigar-op		strand
	with openfile (tsv) as fh:
		for line in fh:
			if line.startswith ('#'):
				continue
			ary = line.strip().split()
			if ary[-2] in ['M','m']:
				k = (ary[2], int (ary[-4]), ary[-1]) #
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				qual[k].append (ord(ary[-5])-33)
				base[k] = ary[-3].upper()
				if (ary[-3] != ary[4]):
					mis[k] += 1
				else:
					mat[k] += 1
			if ary[-2] == 'D':
				k = (ary[2], int(ary[-4]), ary[-1])
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				base[k] = ary[-3].upper() 
				dele[k] = dele.get(k,0) + 1
			if ary[-2] == 'I':
				last_k = aln_mem[-1][1],aln_mem[-1][2],aln_mem[-1][3] # last alignment with match/mismatch/del
				next_k = (ary[2], last_k[1] + 1,last_k[2])
				if last_k[0] != ary[2]:
					sys.stderr.write (line.strip())
				ins_k_up = (ary[0], ary[2], last_k[1],last_k[2])
				ins_k_down = (ary[0], ary[2], last_k[1] + 1,last_k[2])
				if (ins_k_down) not in ins_q:
					ins[next_k] = ins.get(next_k,0) + 1
					ins_q[ins_k_down].append(ord(ary[-5])-33)
				if (ins_k_up) not in ins_q:
					ins[last_k] = ins.get(last_k,0) + 1
					ins_q[ins_k_up].append(ord(ary[-5])-33)
	header = '#Ref,pos,base,cov,mat,mis,ins,del,qual,strand\n'
	#cc6m_2244_T7_ecorv,7,A,1.0,1,0,0,0,15,+
	outh.write(header)
	for k in cov.keys():
		depth = cov.get (k,0) 
		Mis = mis.get (k,0) 
		Mat = mat.get (k,0)
		Del = dele.get (k,0)
		q_lst = qual.get (k,[0])
		try:
			q_lst = ':'.join (map (str, q_lst)) #qual[k]
			num_ins = ins.get (k,0)
			inf = "{},{},{},{},{},{},{},{},{},{}\n".format (k[0], k[1], base[k], depth, Mat, Mis, num_ins, Del, q_lst, k[2])
			outh.write (inf)
		except:
			sys.stderr.write ("file {} {} does not work\n".format (tsv,k))
	return out 

	
def tsv_to_var (tsv):
	'''
	reference base was complemented if aligned on reverse strand 
	'''
	header = "#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del"
	out = '.'.join (tsv.split('.')[:-1]) + '.per.site.var.csv' 
	outh = open (out,'w')
	outh.write('#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del\n')
	mis = defaultdict(int) # mismatches
	mat = defaultdict (int) #matches
	ins = defaultdict(int) # insertions
	dele = defaultdict(int) # deletions
	cov = OrderedDict ()  # coverage
	ins_q = defaultdict(list)
	aln_mem = []  #read, ref, refpos; only store last entry not matching insertion
	pos = defaultdict(list) # reference positions
	base = {} # ref base
	Q = defaultdict(list)
	qual = defaultdict(list)
	basesdict = {'A':'T', 'G':'C','C':'G','T':'A','N':'N'} 
#READ_NAME     FLAG    CHROM   READ_POS  BASE   QUAL  REF_POS REF  OP 	STRAND 
	with openfile (tsv) as fh:
		for line in fh:
			if line.startswith ('#'):
				continue
			ary = line.strip().split()
			if ary[-2] in ['M','m']:
				k = (ary[2], int (ary[-4]), ary[-1]) #
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				qual[k].append (ary[-5])
				Q[k].append(ary[-5])
				base[k] = ary[-3].upper()
				if (ary[-3] != ary[4]):
					mis[k] += 1
				else:
					mat[k] += 1
			if ary[-2] == 'D':
				k = (ary[2], int(ary[-4]), ary[-1])
				cov[k] = cov.get(k,0) + 1
				aln_mem = []
				aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
				base[k] = ary[-3].upper() 
				dele[k] = dele.get(k,0) + 1
			if ary[-2] == 'I':
				last_k = aln_mem[-1][1],aln_mem[-1][2],aln_mem[-1][3] # last alignment with match/mismatch/del
				next_k = (ary[2], last_k[1] + 1,last_k[2])
				if last_k[0] != ary[2]:
					sys.stderr.write (line.strip())
				ins_k_up = (ary[0], ary[2], last_k[1],last_k[2])
				ins_k_down = (ary[0], ary[2], last_k[1] + 1,last_k[2])
				if (ins_k_down) not in ins_q:
					ins[next_k] = ins.get(next_k,0) + 1
					ins_q[ins_k_down].append(ary[-5])
				if (ins_k_up) not in ins_q:
					ins[last_k] = ins.get(last_k,0) + 1
					ins_q[ins_k_up].append(ary[-5])
	#header = '#Ref,pos,base,cov,mis,ins,del,q_sum,strand'
	for k in cov.keys():
		depth = float (cov.get (k,0) )
		Mis = mis.get (k,0) 
		Mat = mat.get (k,0)
		Del = dele.get (k,0)
		q_lst = qual.get (k,[0])
		try:
			num_ins = ins.get (k,0)
			q_mn, q_md, qstd = np.mean (np.array (q_lst).astype(np.float)), np.median (np.array (q_lst).astype(np.float)), np.std (np.array (q_lst).astype(np.float))
			m,i,d = np.array ([Mis,num_ins,Del])/depth
			ref_base = base[k] if k[2] == '+' else basesdict.get (k,'N')
			inf = "{},{},{},{},{},{},{},{},{},{},{}\n".format (k[0], k[1], ref_base, k[2], depth,q_mn,q_md, qstd,m,i,d)
			outh.write (inf)
		except:
			raise 
			sys.stderr.write ("problematic {} in {}\n".format (k,tsv))
	return out
	
def combine_freq (list_of_freq_files):
	'''
	reference base was complemented if reads aligned on reverse strand
	cc6m_2244_T7_ecorv,31,A,4.0,4,0,0,0,4:16:7:3,+
	'''
	mem = defaultdict(lambda: defaultdict(list))
	Var=defaultdict (list)
	Qual = defaultdict (list) 
	ks = OrderedDict()
	outfile = os.path.dirname (list_of_freq_files[0]) if os.path.dirname (list_of_freq_files[0]) else '_'
	outfile=outfile.replace ('_tmp_splitted','') +'.per.site.var.csv' 
	outh = open (outfile,'w')
	basesdict = {'A':'T', 'G':'C','C':'G','T':'A','N':'N'}
	for l in fileinput.input(list_of_freq_files):
		ary = l.strip().split(',')
		k = (ary[0],ary[1],ary[2],ary[-1])
		ks[k] = True
		c,_,m,i,d = map (float,ary[3:8])
		c_m_i_d = np.array ([c,m,i,d])
		q_lst = []
		q_lst = [0] if (len(ary[8])) == 0 else [ float (x) for x in ary[8].split(':')]
		Var[k] = Var.get(k,np.array([0])) + c_m_i_d 
		Qual[k] = Qual.get(k,[]) + q_lst 
	outh.write ('#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del\n')
	basesdict = {'A':'T', 'G':'C','C':'G','T':'A','N':'N'}
	for k in ks:
		cov = Var[k][0]
		q_lst = Qual[k]
		var_freq = Var[k][1:]/cov 
		var_freq = ",".join (var_freq.astype (str))
		k = list (k)
		if k[3] == '-':
			k[2] = basesdict.get(k[2],'N')
		outh.write ("{},{},{},{},{},{}\n".format (','.join(k),cov,'%0.5f'%np.mean(q_lst),'%0.5f'%np.median(q_lst),'%0.5f'%np.std(q_lst), var_freq)) 
	outh.close()
	return (outfile)



def slide_per_site_var_for_unsorted_data (per_site_var,win=5):
	'''
	#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del
	cc6m_2244_T7_ecorv,7,A,+,1.0,15.0,15.0,0.0,0.0,0.0,0.0
	kmer sequences will be reversed if reads aligned on the minus strand 
	bases mapped to reverse strand have alredy been complemented during above processing 
	
	consume a lot of ram for unsorted data 
	'''
	mem = {}
	contents = OrderedDict()
	dist   = int (win)//2 + 1
	fh = open (per_site_var,'r')
	for line in fh:
		if line.startswith ('#'):
			continue
		if re.match ('\s+',line):
			continue
		ary = line.strip().split(',')
		ref,pos,strand = (ary[0],ary[1],ary[3])
		contents[(ref,pos,strand)] = line.rstrip()
	prefix =  per_site_var.replace ('.per_site.var.csv','') # ".".join (per_site_var.split('.')[:-1])
	out_tmp = prefix +'.per.site.var.{}mer.tmp'.format(win)
	outh1= open (out_tmp,'w')
	header = '#Kmer,window,Relative_pos,Ref,Ref_Pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del'
	outh1.write (header+'\n')

	for k in contents.keys():
		ref, pos, strand = k 
		try:
			pos = int(pos)
		except:
			sys.stderr.write("wrong ref pos {}".format(k))
			continue
		POS = []
		LINES = []
		upper = ''
		down = ''
		for i in list (reversed (range(1,dist))):
			POS.append (str(pos-i))
			kk = (ref,str(pos-i),strand)
			if kk in contents:
				base = contents[kk].split(',')[2]
				upper += base
				LINES.append ('-'+str(i)+','+ contents[kk])
			else:
				upper += 'N'
				LINES.append('-'+str(i)+','+'Null')
		LINES.append ('+0'+','+contents[k])
		POS.append (str(pos))
		for j in range (1,dist):
			POS.append (str(pos + j))
			kk = (ref,str(pos+j),strand)
			if kk in contents:
				base = contents[kk].split(',')[2]
				down += base
				LINES.append ('+'+str(j)+','+ contents[kk])
			else:
				down += 'N'
				LINES.append ('+'+str(j)+','+'None')
		positions = '-'.join ([POS[0],POS[-1]])
		for l in LINES:
			kmer = upper+contents[k].split(',')[2] +down
			kmer = kmer if strand == '+'  else  kmer[::-1]
			outh1.write (kmer+',' +positions+','+ l+'\n')
	outh1.close()
	###    sum up slided per site variants from multiple lines into single lines    ############
	'''
	#Kmer,window,Relative_pos,Ref,Ref_Pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del
	TAGGT,1852:1853:1854:1855:1856,-2,cc6m_2459_T7_ecorv,1852,T,+,6490.0,8.62139,8.00000,4.70584,0.1869029275808937,0.3489984591679507,0.02773497688751926
	TAGGT,1852:1853:1854:1855:1856,-1,cc6m_2459_T7_ecorv,1853,A,+,6508.0,6.96834,6.00000,3.65965,0.21204671173939765,0.13844499078057776,0.019668100799016593
	'''	
	mem_window = defaultdict (defaultdict(list).copy)
	
	# !!!!!!!!!!!!!!!!!!!  should rewrite below
	k_pool = []
	f =  open (out_tmp,'r')
	for l in f:
		if l.startswith ('N'):
			continue
		elif l.startswith ('#'):
			continue
		ary = l.strip().split(',')
		if ary[0].endswith ('N'):
			continue
		if (len(ary) < 14):
			continue
		k = ''
		try:
			ks = (ary[0], ary[1],ary[3],ary[6]) # #Kmer,window,Ref,strand
			k = ",".join (ks)
			k_pool.append (k)
		except:
			pass
		try:
			mem_window[k]['q'].append(ary[8])
		except:
			sys.stderr.write ("problematic line for q: " + l.strip())
		try:
			mem_window[k]['m'].append(ary[11])
		except:
			sys.stderr.write ("problematic line for m" + l.strip())
		try:
			mem_window[k]['i'].append(ary[12])
		except:
			sys.stderr.write ("problematic line for i:" + l.strip())
		try:
			mem_window[k]['d'].append(ary[13])
		except:
			sys.stderr.write ("problematic line for d" + l.strip())
		try:
			mem_window[k]['cov'].append(ary[7])
		except:
			sys.stderr.write ('problematic line for cov:' + l.strip())
	f.close()
	out2 = prefix + '.per.site.{}mer.csv'.format(win)
	outh2 = open (out2,'w')
	q_in_head = ",".join (["q{}".format(i) for i in range(1,win+1)])
	mis_in_head = ",".join (["mis{}".format(i) for i in range(1,win+1)])
	ins_in_head = ",".join (["ins{}".format(i) for i in range(1,win+1)])
	del_in_head = ",".join (["del{}".format(i) for i in range(1,win+1)])
	outh2.write ('#Kmer,Window,Ref,Strand,Coverage,{},{},{},{}\n'.format(q_in_head, mis_in_head, ins_in_head, del_in_head))
	for k in set (k_pool):
			Qs   = ",".join (mem_window[k]['q'])
			Mis  = ",".join (mem_window[k]['m'])
			Ins  = ",".join (mem_window[k]['i'])
			Del  = ",".join (mem_window[k]['d'])
			Cov = ":".join (mem_window[k]['cov'])
			outh2.write (",".join ([k,Cov,Qs,Mis,Ins,Del])+'\n')
	outh2.close()
	os.remove (out_tmp)
	return (out2)
	
def print_last_consecutive_lines (lines, outfh):
	contents = OrderedDict()
	for line in lines:
		ary = line.strip().split(',')
		ref,pos,strand = (ary[0], ary[1], ary[3])
		contents[(ref,pos,strand)] = line.rstrip()
		
	win = len (lines)
	middle = lines [win//2].rstrip().split(',')
	window = str(int(middle[1]) - win//2)+'-'+str(int(middle[1]) + win//2 )
	kmer = ''
	
	consecutive_lines = []
	ref,pos,base,strand = middle[:4]
	for i in reversed (list (range(1, win//2+1))):
		k = (ref,str(int(pos) -i), strand)
		relative_pos = '-'+str(i)
		if k in contents:
			kmer = kmer + contents[k].split(',')[2]
			consecutive_lines.append (window + ',' + relative_pos +','+contents[k])
		else:
			kmer = kmer + 'N'
			consecutive_lines.append (window +','+relative_pos+','+","+",".join (['NA']*12))
			
	consecutive_lines.append (window +',' + '+0' +','+",".join (middle))
	kmer = kmer + middle[2]
	for i in range(1,win//2+1):
		k = (ref,str(int(pos) + i),strand)
		relative_pos = '+'+str(i)
		if k in contents:
			kmer = kmer + contents[k].split(',')[2]
			consecutive_lines.append (window + ',' + relative_pos+','+contents[k])
		else:
			kmer = kmer + 'N'
			consecutive_lines.append (window +',' + relative_pos+','+",".join (['NA']*12))
	for l in consecutive_lines:
		print (kmer+','+l, file=outfh)
		

def slide_per_site_var (per_site_var,win=5):
	'''
	#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del

	cc6m_2244_T7_ecorv,7,A,+,1.0,15.0,15.0,0.0,0.0,0.0,0.0
	kmer sequences will be reversed if reads aligned on the minus strand 
	bases mapped to reverse strand have alredy been complemented during above processing 
	'''
	#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del
	prefix =  re.sub(r'.per.site.\S+','',per_site_var)# , .replace ('.per.site.csv','') # ".".join (per_site_var.split('.')[:-1])
	out_tmp = prefix +'.per_site_var.{}mer.tmp'.format(win)
	if os.path.exists (out_tmp):
		os.remove (out_tmp)
	outfh = open (out_tmp,'w')
	
	fh = open (per_site_var, 'rb' )
	eof = fh.seek (-1,2)
	fh.seek(0,0)
	head = fh.readline () 
	lines = [] 

	for _ in range (win):
		l = fh.readline().decode('utf-8').rstrip()
		if l:
			lines.append (l)
	if len  (lines)	 < win:
		print ('not enough sites to be slided',file=sys.stderr)
	
	contents = OrderedDict()

	for line in lines:
		ary = line.strip().split(',')
		ref,pos,strand = (ary[0], ary[1], ary[3])
		contents[(ref,pos,strand)] = line.rstrip()
	
	while (fh.tell() <= eof):
		middle = lines [win//2].split(',')
		window = str(int(middle[1]) - win//2)+'-'+str(int(middle[1]) + win//2 )
		consecutive_lines = []
		kmer=''
		ref,pos,base,strand = middle[:4]
		k_to_del = (ref,str(int(pos)-win),strand)
		for i in reversed (list (range(1, win//2+1))):
			k = (ref,str(int(pos) -i),strand)
			relative_pos = '-'+str(i)
			if k in contents:
				kmer = kmer +contents[k].split(',')[2]
				consecutive_lines.append (window+','+relative_pos +','+contents[k])
			else:
				consecutive_lines.append (window+','+relative_pos+','+ "," . join ([ref, str(int(pos) -i), 'N', strand, '0', 'NaN,NaN,NaN,NaN,NaN,NaN']))
				kmer = kmer+'N'
		consecutive_lines.append (window+',+0' +','+",".join (middle))
		kmer = kmer +middle[2]
		for i in range(1,win//2+1):
			k = (ref,str(int(pos)+i),strand)
			relative_pos = '+'+str(i)
			if k in contents:
				kmer = kmer +contents[k].split(',')[2]
				consecutive_lines.append (window+','+relative_pos+','+contents[k])
			else:
				kmer = kmer+'N'
				consecutive_lines.append (window+','+relative_pos+','+ "," . join ([ref, str(int(pos) +i), 'N', strand, '0', 'NaN,NaN,NaN,NaN,NaN,NaN']))
				#consecutive_lines.append (window+','+relative_pos+','+ "," . join (['NaN']*11))
				
		for l in consecutive_lines:
			print (kmer+','+l,file = outfh)
		keys = list(contents.keys())
		del consecutive_lines
		if k_to_del in contents:
			del contents[k_to_del]
		lines = lines[1:]
		new_line = fh.readline().decode('utf-8').rstrip()
		lines.append (new_line)
		ref,pos,base,strand = new_line.split(',')[:4]
		contents[(ref,pos,strand)] = new_line
	print_last_consecutive_lines (lines, outfh)
	outfh.close()
	
	#out2 = prefix + '.per_site.{}mer.csv'.format(win)
	out2 = prefix + '.per.site.{}mer.csv'.format(win)
	outh2 = open (out2,'w')
	q_in_head = ",".join (["q{}".format(i) for i in range(1,win+1)])
	mis_in_head = ",".join (["mis{}".format(i) for i in range(1,win+1)])
	ins_in_head = ",".join (["ins{}".format(i) for i in range(1,win+1)])
	del_in_head = ",".join (["del{}".format(i) for i in range(1,win+1)])
	outh2.write ('#Kmer,Window,Ref,Strand,Coverage,{},{},{},{}\n'.format(q_in_head, mis_in_head, ins_in_head, del_in_head))
	
	tmpfh = open (out_tmp,'r')
	cov, q, mis, ins, dele = [], [], [], [], []
	firstline = tmpfh.readline().rstrip().split(',')
	current_win = (firstline[0], firstline[1], firstline[3], firstline[6])
	lines = []
	lines.append (firstline)
	ary = []
	for l in tmpfh:
		ary = l.rstrip().split(',')
		try:
			window = (ary[0], ary[1], ary[3], ary[6])		
		except:
			print (l.rstrip())
		if window != current_win:
			for ele in lines:
				q.append (ele[8])
				mis.append (ele[11])
				ins.append (ele[12])
				dele.append (ele[13])
				cov.append(ele[7])
			Qs = ",".join (q)
			Mis = ",".join (mis)
			Ins = ",".join (ins)
			Del = ",".join(dele)
			Cov = ":".join (cov)
			print (",".join (current_win), Cov, Qs, Mis, Ins, Del, sep=",", file= outh2)
			cov, q, mis, ins, dele = [], [], [], [], []
			current_win = window
			lines = []
		lines.append (ary)
	# last 5 lines
	cov, q, mis, ins, dele  = [], [], [], [], []
	for ele in lines:
		q.append (ele[8])
		mis.append (ele[11])
		ins.append (ele[12])
		dele.append (ele[13])
		cov.append (ele[7])
	Qs = ",".join (q)
	Mis = ",".join (mis)
	Ins = ",".join (ins)
	Del = ",".join(dele)
	Cov = ":".join (cov)
	print (",".join (window), Cov, Qs, Mis, Ins, Del, sep=",", file= outh2)
	
	tmpfh.close()
	outh2.close()
	os.remove (out_tmp)
	return (out2)
	
	
def per_read_var (tsv):
	''' single thread'''
	rdnames = []
	qualities = {}
	dels = defaultdict(int)
	mis = {}
	ins = {}
	match = {}
	ref_pos = {}
	k = ''
	next_k = ''
	
	if tsv.endswith (".gz"):
		fh = gzip.open (tsv)
	else:
		fh = open (tsv)
	for line in fh:
		if re.match ('\s+',line):
			continue
		if re.match('#',line):
			continue
		if re.match (':',line):
			continue
		ary = line.rstrip().split()
		try:
			if not re.match (r'[MID]',ary[8]):
				continue
		except:
			print ('problematic line:', line)
		if len (ary) != 10:
			continue
		if ary[6].startswith('-'):
			continue
		if re.match (r'[HS]',ary[8]):
			continue
	#READ_NAME FLAG CHROM READ_POS BASE QUAL REF_POS REF OP
		if ary[6] != '.':
			ary[6] = str(int(ary[6])) # ref pos is 1-based
		ref = ary[2]
		ref_pos = ary[6]
		ref_base = ary[7]
		rd = ary[0]
		rd_pos = ''
		if ary[8] != 'D':
			rd_pos = str (int (ary[3]) + 1)#turn read_pos into 1-based
			rd_base = ary[4]
			strand = ary[-1]
			k = ','.join ([ref,ref_pos,strand,rd,rd_pos]) #include reference to account for multi-mappings
			qualities[k] = ord(ary[5]) - 33
			dels[k] = dels.get(k,0) + 0
			rdnames.append(k)
			if ary[8].upper() == 'M' and ary[4] == ary[7]: # or re.match ('I', ary[-1].upper()):
				mis[k] = '0'
				ins[k] = '0'
			elif ary[8].upper() == 'M' and ary[4] != ary[7]: # or re.match ('I', ary[-1].upper()):
				mis[k] = '1'
				ins[k] = '0'
			elif ary[8].upper() == 'I':
				ins[k] = '1'
				mis[k] = '0'
		elif ary[8].upper() == 'D':
			dels[k]= dels.get(k,0) + 1
			mis[k] = '0'
			ins[k] = '0'
		else:
			continue
	prefix = re.sub(r'.tsv$', '', tsv)
	del_tmp = prefix + '.per_read_var.tmp.csv'
	tmp_fh = open (del_tmp,'w')
	if tsv.endswith (".gz"):
		fh = gzip.open (tsv)
	else:
		fh = open (tsv)
	for l in fh:
		if l.startswith ('#'):
			continue
		if re.match ('\s+',l):
			continue
		if re.match('#',l):
			continue
		if re.match (':',l):
			continue
		ary = l.rstrip().split()
		if len (ary) != 10:
			continue
		if ary[6].startswith('-') :
			continue
		ref = ary[2]
		ref_pos = ary[6]
		ref_base = ary[7]
		rd = ary[0]
		rd_pos = '' #ary[3]
		strand = ary[-1]
		if not re.match (r'[DMI]',ary[8]):
			continue
		if all ([ary[8] != 'D' , ary[8] != 'S' , ary[8] != 'H']):
			rd_pos = str (int (ary[3]) + 1)
			rd_base = ary[4]
			k = ','.join ([ref,ref_pos, strand,rd,rd_pos])
			inf = k.split(',')
			tmp_fh.write (','.join ([inf[0],inf[1],ref_base,inf[2],inf[3],inf[4],rd_base,str 
			(qualities[k]),mis[k],ins[k],str(dels[k])])+'\n')
	tmp_fh.close()
	
	adjusted_file = prefix+'_per_rd_var.csv'
	adj_fh = open (adjusted_file,'w')
	adj_fh.write ( ",".join (
		["#REF",'REF_POS','REF_BASE','STRAND','READ_NAME','READ_POSITION','READ_BASE','BASE_QUALITY','MISMATCH','INSERTION','DELETION']) + '\n')
	last_del = defaultdict(float)
	line_num = 0
	with open (del_tmp,'r') as FH:
		for l in FH:
			ary = l.strip().split(',')
			half = float(ary[-1]) / 2
			ary[-1] = last_del.get (line_num,0) + half
			line_num += 1
			last_del[line_num] = half
			adj_fh.write ( ','.join(map (str, ary)) + '\n')
	adj_fh.close()
#	os.remove (del_tmp)
	return adjusted_file

def split_tsv_for_per_read_var (tsv, q, number_threads):
	''' put single read into queue collector '''
	single_read_tsv_lines = []  # store single read rather than multiplt read 
	last_seen = ''
	read_cnt = 0
	'''
	if isinstance (tsv,str):  #intput is tsv text file 
		fh = openfile (tsv)  
	elif isinstance (tsv,subprocess.Popen): #return from subproces.Popen 
		fh = tsv.stdout
	elif hasattr(tsv, 'read'):
		fh = tsv 
	'''
	try:
		if True:
			for l in tsv:
				if isinstance (l, bytes):
					l = l.decode ('utf-8')
				if l and l.startswith ('#'):
					continue
				rd = l.split()[0]
				if rd != last_seen:  #, single_read_lines, read_cnt>=1]):
					last_seen = rd
					if read_cnt > 0:
						q.put( (read_cnt, single_read_tsv_lines))
						single_read_tsv_lines = []
						read_cnt += 1
					else:  # first read 
						read_cnt += 1 
				single_read_tsv_lines.append (l.strip())
		q.put((read_cnt, single_read_tsv_lines))
		sys.stderr.write("totally processed {} reads\n".format(read_cnt))
	except:
		raise 
		sys.stderr.write("split tsv file on reads failed\n")
	finally:
		for _ in range(number_threads ):
			q.put(None)

def per_read_var_multiprocessing (inq, number_threads, outputfile):
	outfh = open (outputfile, 'a')
	for _, single_read_lines in iter (inq.get, None):
		copy_single_read_tsv_lines = []
		qualities = {}
		dels = defaultdict(int)
		mis = {}
		ins = {}
		match = {}
		ref_pos = {}
		k = ''
		next_k = ''
		for line in single_read_lines:
			copy_single_read_tsv_lines.append (line)
			if re.match ('\s+',line):
				continue
			if re.match('#',line):
				continue
			if re.match (':',line):
				continue
			ary = line.rstrip().split()
			try:
				if not re.match (r'[MID]',ary[8]):
					continue
			except:
				print ('problematic line:', line)
			if len (ary) != 10:
				continue
			if ary[6].startswith('-'):
				continue
			if re.match (r'[HS]',ary[8]):
				continue
			if ary[6] != '.':
				ary[6] = str(int(ary[6])) # ref pos is 1-based
			ref = ary[2]
			ref_pos = ary[6]
			ref_base = ary[7]
			rd = ary[0]
			rd_pos = ''
			if ary[8] != 'D':
				rd_pos = str (int (ary[3]) + 1)#turn read_pos into 1-based
				rd_base = ary[4]
				strand = ary[9]
				k = ','.join ([ref,ref_pos,strand,rd,rd_pos]) #include reference to account for multi-mappings
				qualities[k] = ord(ary[5]) - 33
				dels[k] = dels.get(k,0) + 0
				#rdnames.append(k)
				if ary[8].upper() == 'M' and ary[4] == ary[7]: 
					mis[k] = '0'
					ins[k] = '0'
				elif ary[8].upper() == 'M' and ary[4] != ary[7]:
					mis[k] = '1'
					ins[k] = '0'
				elif ary[8].upper() == 'I':
					ins[k] = '1'
					mis[k] = '0'
			elif ary[8].upper() == 'D':
				dels[k]= dels.get(k,0) + 1
				mis[k] = '0'
				ins[k] = '0'
			else:
				continue
		single_read_var_info = []
		for l in copy_single_read_tsv_lines:
			if l.startswith ('#'):
				continue
			if re.match ('\s+',l):
				continue
			if re.match('#',l):
				continue
			if re.match (':',l):
				continue
			ary = l.rstrip().split()
			if len (ary) != 10:
				continue
			if ary[6].startswith('-') :
				continue
			ref = ary[2]
			ref_pos = ary[6]
			ref_base = ary[7]
			rd = ary[0]
			rd_pos = '' #ary[3]
			strand = ary[-1]
			if not re.match (r'[DMI]',ary[8]):
				continue
			if all ([ary[8] != 'D' , ary[8] != 'S' , ary[8] != 'H']):
				rd_pos = str (int (ary[3]) + 1)
				rd_base = ary[4]
				k = ','.join ([ref,ref_pos, strand,rd,rd_pos])
				inf = k.split(',')
				single_read_var_info.append ([inf[0],inf[1],ref_base,inf[2],inf[3],inf[4],rd_base,str 
				(qualities[k]),mis[k],ins[k],str(dels[k])])
		last_del = defaultdict(float)
		line_num = 0
		for ary in single_read_var_info:
			half = float(ary[-1]) / 2
			ary[-1] = last_del.get (line_num,0) + half
			line_num += 1
			last_del[line_num] = half
			outfh.write ( ','.join(map (str, ary)) + '\n')
			
def split_reads_for_per_read_var_sliding (per_read_var,q,number_threads):
	#REF,REF_POS,REF_BASE,STRAND,READ_NAME,READ_POSITION,READ_BASE,BASE_QUALITY,MISMATCH,INSERTION,DELETION
	last_seen = ''
	rd = '' 
	read_cnt = 0 
	single_read_var_lines = []
	try:
		with openfile (per_read_var) as fh:
			for l in fh:
				if l and l.startswith ('#'):
					continue
				rd = l.split(',')[4]
				if rd != last_seen:  #, single_read_lines, read_cnt>=1]):
					last_seen = rd
					if read_cnt > 0:
						q.put( (read_cnt, single_read_var_lines))
						single_read_var_lines = []
						read_cnt += 1
					else:  # first read 
						read_cnt += 1 
				single_read_var_lines.append (l.strip())
		q.put((read_cnt, single_read_var_lines))
		sys.stderr.write("totally processed {} reads\n".format(read_cnt))
	except:
		raise 
		sys.stderr.write("split tsv file on reads failed\n")
	finally:
		for _ in range(number_threads ):
			q.put(None)
		
def slide_per_read_var_multiprocessing (inq, outfile):
	#REF,REF_POS,REF_BASE,STRAND,READ_NAME,READ_POSITION,READ_BASE,BASE_QUALITY,MISMATCH,INSERTION,DELONI
	for idx, read_var_lines in iter(inq.get,None):
		relative_pos = ['-2','-1','0','+1','+2']
		rd_rf_mem_di = [] #defaultdict(list)
		ref_strand = []
		rd_bases = [] #defaultdict(str)
		rf_bases = [] #defaultdict(str)
		rd_pos = [] #defaultdict(str)
		rf_pos = [] #defaultdict(str)
		qual = [] #defaultdict(str)
		mis = [] #defaultdict(str)
		ins = [] #defaultdict(str)
		delt = [] #defaultdict(str)
		prefix = outfile  +'.'+str(idx)
		row_counter = 1
		slided_file = prefix+'.5mer.tmp.csv'
		slided_fh = open (slided_file,'w')
		slided_header = '#Relative_Pos,window,ReadKmer,RefKmer,Ref,RefPos,RefBase,Strand,Read,ReadPos,ReadBase,Q,M,I,D'
		slided_fh.write (slided_header+'\n')
		#REF,REF_POS,REF_BASE,STRAND,READ_NAME,READ_POSITION,READ_BASE,BASE_QUALITY,MISMATCH,INSERTION,DELONI
		# first line always starts with a read containing perfect match or variant other than deletion;
		# slided var
		firstLine = read_var_lines[0]
		tmp = firstLine.strip().split(',')
		rd = tmp[4]
		ref = tmp[0]
		strand = tmp[3]
		comb = rd+' '+ref+' ' + strand
		ref_strand.append (strand)
		rd_rf_mem_di.append (comb)
		rd_bases.append (tmp[6])
		rf_bases.append (tmp[2])
		rd_pos.append (tmp[5])
		rf_pos.append (tmp[1])
		qual.append(tmp[7])
		mis.append(tmp[8])
		ins.append(tmp[9])
		delt.append(tmp[10])
		
		for line in read_var_lines[1:]:
			if line.startswith('#'):
				continue
			tmp = line.strip().split(',')
			if len (tmp) < 10:
				continue
			rd = tmp[4]
			ref = tmp[0]
			comb = rd+' '+ref + ' ' + tmp[3]
			if comb in rd_rf_mem_di:
				rd_rf_mem_di.append (comb)
				rd_bases.append (tmp[6])
				rf_bases.append (tmp[2])
				rd_pos.append (tmp[5])
				rf_pos.append (tmp[1])
				qual.append(tmp[7])
				mis.append(tmp[8])
				ins.append(tmp[9])
				delt.append(tmp[10])
			else:
				for i in window(range(len(rd_rf_mem_di)), 5):
					last = None
					rd_win_pos = ''
					if i[-1] is None: # when a read have less than window size number of positions mapped to a reference
						last = i[0] + len (i) + 1
						rd_win_pos = ':'.join(rd_pos[i[0]:last])
					else:
						last = i[-1] + 1
						rd_win_pos = ':'.join(rd_pos[i[0]:last])
				#sys.stderr.write ('exception',i, rd_rf_mem_di[i], rd_pos[i],rf_pos[i])
						#sys.stderr.write ('exception',i, rd_rf_mem_di[i], rd_pos[i],rf_pos[i])
					RD, REF,STRAND = rd_rf_mem_di[-1].split()
					rd_kmer = ''.join (rd_bases[i[0]:last])
					rf_kmer = ''.join (rf_bases[i[0]:last])
					int_ele = 0
					for k in (i):
						if isinstance (k,int):
					#sys.stderr.write (str (len (rd_rf_mem_di)) + ' rf_rd: ' + rd_rf_mem_di[0] + ' ' + 
					#rd_rf_mem_di[-1] + ' '+ ':'.join (map (str, rd_pos)) + ' ' + ':'.join (map (str, 
					#i))+'\n' )
							slided_fh.write (",".join ([relative_pos[int_ele], rd_win_pos, rd_kmer, rf_kmer, REF,  rf_pos[k],rf_bases[k], STRAND,
										RD,rd_pos[k],rd_bases[k],
										qual[k],mis[k],ins[k],delt[k]]) + '\n')
							int_ele += 1
						elif k is None: # less than 5 sites of a read mapped
							slided_fh.write (",".join ([relative_pos[int_ele], rd_win_pos, rd_kmer, rf_kmer, REF, 'Nan,Nan,Nan',
										RD,'Nan,Nan,Nan,Nan,Nan,Nan'])+'\n')
							int_ele += 1
				rd_rf_mem_di = [] #defaultdict(list)
				rd_bases = [] #defaultdict(str)
				rf_bases = [] #defaultdict(str)
				rd_pos = [] #defaultdict(str)
				rf_pos = [] #defaultdict(str)
				qual = [] #defaultdict(str)
				mis = [] #defaultdict(str)
				ins = [] #defaultdict(str)
				delt = [] #defaultdict(str)
				rd_rf_mem_di.append (comb)
				rd_bases.append (tmp[6])
				rf_bases.append (tmp[2])
				rd_pos.append (tmp[5])
				rf_pos.append (tmp[1])
				qual.append(tmp[7])
				mis.append(tmp[8])
				ins.append(tmp[9])
				delt.append(tmp[10])
		
		for i in window(range(len(rd_rf_mem_di)),5):
			last = None
			rd_win_pos = ''
			if i[-1] is None:
				last = i[0] + len(i) + 1
				rd_win_pos = ':'.join(rd_pos[i[0]:last])
			else:
				last = i[-1] + 1
				rd_win_pos = ':'.join(rd_pos[i[0]:last])
			RD, REF,STRAND = rd_rf_mem_di[-1].split()
			rd_kmer = ''.join (rd_bases[i[0]:last])
			rf_kmer = ''.join (rf_bases[i[0]:last])
			int_ele = 0
			for k in (i):
				if isinstance (k, int):
					slided_fh.write (",".join ([relative_pos[int_ele],rd_win_pos,rd_kmer, rf_kmer, REF, rf_pos[k], rf_bases[k], STRAND,
										RD,rd_pos[k],rd_bases[k],qual[k],mis[k],ins[k],delt[k]]) + '\n')
					int_ele += 1
				elif k is None:
					slided_fh.write (",".join ([relative_pos[int_ele],rd_win_pos,rd_kmer, rf_kmer, REF, 'Nan,Nan,Nan',
										RD,'Nan,Nan,Nan,Nan,Nan,Nan']) + '\n')
					int_ele += 1
		
		slided_fh.close()
		'''
		-2,80:81:82:83:84,GTTTT,GTTTT+C+G,cc6m_2244_T7_ecorv,114,G,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,80,G,6,0,0,0
		...
		+2
		'''
		##########################################################################
		##############  sum up onebase oneline into onekmer one line format ######
		##########################################################################
		
		sum_out = open (outfile,'a')
		
		'''
		input
		#Relative_Pos,window,ReadKmer,RefKmer,Ref,RefPos,RefBase,Strand,Read,ReadPos,ReadBase,Q,M,I,D
		-2,2:3:4:5:6,TCGAT,TCGAT,cc6m_2244_T7_ecorv,46,T,94e2685c-c9e3-42b5-9a36-2c9ed081c1bc,2,T,5,0,0,0
		...
		+2
		'''
		from collections import defaultdict
		mem_window = defaultdict (defaultdict(list).copy)
		f = openfile (slided_file)
		header = f.readline()
		#new_header = ",".join (header.strip().split(',')[2:])
		#sum_out.write ("#Read,Read_Window,ReadKmer,Ref,RefKmer,Strand,Ref_Window,q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,ins1,ins2,ins3,ins4,ins5,del1,del2,del3,del4,del5\n")
		ary = f.readline().strip().split(',')
		k = ','.join(ary[1:5]) + ','+ary[7]+','+ary[8]
		mem_window[k]['q'].append(ary[-4])
		mem_window[k]['m'].append(ary[-3])
		mem_window[k]['i'].append(ary[-2])
		mem_window[k]['d'].append(ary[-1])
		mem_window[k]['refpos'].append(ary[5])
		
		for l in f:
			ary = l.strip().split(',')
			k = ','.join(ary[1:5]) + ','+ary[7]+','+ary[8]
			if k in mem_window:
				mem_window[k]['q'].append(ary[-4])
				mem_window[k]['m'].append(ary[-3])
				mem_window[k]['i'].append(ary[-2])
				mem_window[k]['d'].append(ary[-1])
				mem_window[k]['refpos'].append(ary[5])
			else:
				for kk in mem_window:
					lst = kk.split(',')
					h = ",".join ([lst[5],lst[0],lst[1],lst[3],lst[2],lst[4]])
					sum_out.write (','.join ([h ,':'.join (mem_window[kk]['refpos']), ','.join (mem_window[kk]['q']),
						','.join (mem_window[kk]['m']), ','.join (mem_window[kk]['i']), 
						','.join (	mem_window[kk]['d'])]) +'\n')
				del mem_window[kk]
				mem_window[k]['q'].append(ary[-4])
				mem_window[k]['m'].append(ary[-3])
				mem_window[k]['i'].append(ary[-2])
				mem_window[k]['d'].append(ary[-1])
				mem_window[k]['refpos'].append(ary[5])
		for kk in mem_window.keys():
			lst = kk.split(',')
			h = ",".join ([lst[4],lst[0],lst[1],lst[3],lst[2],lst[4]])
			sum_out.write (','.join ([h ,':'.join (mem_window[kk]['refpos']), ','.join (mem_window[kk]['q']), 
								','.join (mem_window[kk]['m']), ','.join (mem_window[kk]['i']), 
								','.join (mem_window[kk]['d'])]) + '\n')
		f.close()
		os.remove (slided_file)
		
		
		
def _split_eventalign_tbl_on_read (fh , in_q, tmp_dir, num_cpus=3, num_reads_in_each_chunk=4000):
	'''
	split nanopolish evntalign reults on reads\
	and 	collapse on kmer the current intensity and dwell time details.
	...
	Parameters
	~~~~~~~~~
	eventalign_fn: nanopolish event align results, generated with command:
	 			   [nanopolish eventalign  --reads  ${fastq} --bam ${bam} --genome ${ref}\
				    --scale-events  -t 6 -n --samples --signal-index]
	num_reads_in_each_chunk: will store as many as this nuber of reads in each chunk
	in_q: multiprocess Queue to store chunks
	num_cpus: how_many_threads will be used to further process the small chunks after splitting

	Return
	---------
	in_q is return in the end to be processed by downstream functions

	Nanopolish	eventalign results format
	--------------------------------------
	contig  position        reference_kmer  read_name\
	strand  event_index     event_level_mean        event_stdv\
	event_length model_kmer      model_mean      model_stdv\
	standardized_level      start_idx       end_idx samples
	'''
	
	try:
		small_chunk = []
		idx = 0
		reads_cnt = 0
		new_start = ''
		header = ''
		current_read = ''
		smallfile = f"{tmp_dir}/{idx}.chunk"
		outfh = open (smallfile, 'w')
		new_start_read = ''
				
		for line in fh:
			ary = line.rstrip().split()
			if line and line.startswith ('contig'):
				header = "\t".join ([ary[0], ary[1],ary[2],ary[4],ary[8],ary[9],ary[15]])
				print (header, file = outfh)
				continue
			rd = ary[3]
			line = "\t".join ([ary[0], ary[1],ary[2],ary[4],ary[8],ary[9],ary[15]]) #columns: 1-3,5,9-10,16)
			if current_read != rd:
				reads_cnt += 1
				current_read = rd
			if (reads_cnt >= num_reads_in_each_chunk and reads_cnt % num_reads_in_each_chunk == 1):
				if new_start_read != current_read:
					idx += 1
					in_q.put (smallfile)
					outfh.close()
					smallfile = f"{tmp_dir}/{idx}.chunk"
					outfh = open (smallfile,'w')
					print (header, file = outfh)
					new_start_read = current_read
			print (line, file=outfh)
		in_q.put (smallfile)
		outfh.close()
	except Exception:
		print ('splittting nanopolish eventalign results failed',file=sys.stderr)
		raise
	finally:
		print ('finish splitting',file = sys.stderr)
		for i in range (num_cpus):
			in_q.put (None)

def load_single_file (file):
	'''
	load single nanopolish event align file to dataframe 
	'''
	columns ="contig  position reference_kmer  strand model_kmer event_length samples".split()
	if file.endswith ('.gz'):
		return (dd.read_csv (file, sep="\s+",usecols=columns, compression='gzip'))
	else:
		return (dd.read_csv (file, sep="\s+",usecols=columns))




if __name__ == '__main__':
	main ()


