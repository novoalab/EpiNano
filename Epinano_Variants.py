#!/usr/bin/env python 
# -*- coding: utf-8 -*- 

import sys,os,re,io
import shutil, fileinput 
import glob, itertools 
import subprocess 
import argparse
import multiprocessing as mp
from multiprocessing import Process, Manager 
from sys import __stdout__
import dask.dataframe as dd 
from collections import defaultdict 
from collections import OrderedDict 
import numpy as np 

#~~~~~~~~~~~~~~~~~~~~ private function ~~~~~~~~
# func1 subprocess call linux cmmands 

def openfile(f):
	if f.endswith ('.gz'):
		fh = gzip.open (f,'rt')
	elif f.endswith ('bz') or f.endswith ('bz2'):
		fh = bz2.open(f,'rt')
	else:
		fh = open(f,'rt')
	return fh

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



def split_tsv_for_per_site_var_freq(tsv, folder, q, number_threads, num_reads_per_chunk=4000):
	head = next(tsv)
	firstline = next (tsv)
	current_rd = firstline.split()[0]
	rd_cnt = 1
	idx = 0
	out_fn = "{}/CHUNK_{}.txt".format(folder, idx)
	out_fh = open (out_fn,'w')
	print (firstline.rstrip(), file=out_fh)
	try:
		for line in tsv:
			rd = line.split()[0]
			if current_rd != rd:
				rd_cnt += 1
				current_rd = rd
				if ((rd_cnt-1) % num_reads_per_chunk == 0 and rd_cnt >= num_reads_per_chunk):
					out_fh.close()
					q.put ((idx, out_fn))
					idx += 1
					out_fn = "{}/CHUNK_{}.txt".format(folder,idx)
					out_fh = open (out_fn, 'w')
			print (line.rstrip(), file = out_fh)
		out_fh.close()
		q.put ((idx, out_fn))
	except:
		raise
		sys.stderr.write("split tsv file on reads failed\n")
	finally:
		for _ in range(number_threads):
			q.put(None)

def file_exist (file):
	return os.path.exists (file)
	
def _rm (file):
	os.remove (file)
	
def stdin_stdout_gen (stdin_stdout):
	'''
	generator for subprocess popen stdout 
	'''
	for l in stdin_stdout:
		if isinstance (l,bytes):
			yield (l.decode('utf-8'))
		else:
			yield l 

def java_bam_to_tsv (bam_file,  reference_file, sam2tsv, type):
	'''
	type: reference types,i.e., trans or genome 
	'''
	awk_forward_strand = """ awk '{if (/^#/) print $0"\tSTARAND"; else print $0"\t+"}' """
	awk_reverse_strand = """ awk '{if (/^#/) print $0"\tSTARAND"; else print $0"\t-"}' """
	cmds = []

	if type.lower().startswith ("t"):	
		cmd =  f"samtools view -h -F 3860 {bam_file} | java -jar  {sam2tsv} -r {reference_file} "\
			f" | {awk_forward_strand}"		
		cmds = [cmd]
	else:
		cmd1 = (f"samtools view -h -F 3860 {bam_file} | java -jar  {sam2tsv} -r {reference_file} "
					f"| {awk_forward_strand} ")
		cmd2 = (f"samtools view -h -f 16 -F 3844 {bam_file} | java -jar  {sam2tsv} -r {reference_file} "
					f" | {awk_reverse_strand}")	
		cmds = [cmd1,cmd2]
	return cmds 
# data frame 
def tsv_to_freq_multiprocessing_with_manager (tsv_reads_chunk_q, out_dir):
	'''
		mutliprocessing
		produced with sam2tsv.jar with strand information added
		read read-flags reference       read-pos        read-base       read-qual       ref-pos ref-base                cigar-op                strand
		a3194184-d809-42dc-9fa1-dfb497d2ed6a    0       cc6m_2244_T7_ecorv      0       C       #       438     G       S       +
	'''
	for idx, tsv_small_chunk in iter (tsv_reads_chunk_q.get, None):
		filename = "{}/small_{}.freq".format(out_dir, idx)
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
		read_bases = defaultdict (dict)
                #READ_NAME     FLAG    CHROM   READ_POS  BASE   QUAL  REF_POS REF  OP   STRAND
                #read read-flags        reference       read-pos        read-base       read-qual       ref-pos ref-base                cigar-op                strand
		#print ("tsv to freq for chunking",idx)
		with open (tsv_small_chunk, 'r') as fh:
			for line in fh:
				if line.startswith ('#'):
					continue
				ary = line.rstrip().split()
				if ary[-2] in ['M','m']:
					ary[4] = ary[4].upper() # in case soft masked
					ary[7] = ary[7].upper()
					k = (ary[2], int (ary[-4]), ary[-1]) #
					cov[k] = cov.get(k,0) + 1
					aln_mem = []
					aln_mem.append((ary[0],ary[2],int(ary[-4]), ary[-1]))
					qual[k].append (ord(ary[-5])-33)
					base[k] = ary[-3].upper()
					read_bases[k][ary[4]] = read_bases[k].get(ary[4], 0) + 1
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
						pass
					ins_k_up = (ary[0], ary[2], last_k[1],last_k[2])
					ins_k_down = (ary[0], ary[2], last_k[1] + 1,last_k[2])
					if (ins_k_down) not in ins_q:
						ins[next_k] = ins.get(next_k,0) + 1
						ins_q[ins_k_down].append(ord(ary[-5])-33)
					if (ins_k_up) not in ins_q:
						ins[last_k] = ins.get(last_k,0) + 1
						ins_q[ins_k_up].append(ord(ary[-5])-33)
		header = '#Ref,pos,base,cov,mat,mis,ins,del,qual,strand,bases\n'
		outh.write(header)
		os.remove (tsv_small_chunk)	
		for k in cov.keys():
			depth = cov.get (k,0)
			Mis = mis.get (k,0)
			Mat = mat.get (k,0)
			Del = dele.get (k,0)
			q_lst = qual.get (k,[0])
			try:
				q_lst = ':'.join (map (str, q_lst))+':'  # dataframe sum
				num_ins = ins.get (k,0)
				bases_counts = "0:0:0:0:"
				if k in read_bases:
					bases_counts = ":".join ([str(read_bases[k].get(l,0)) for l in 'ACGT'])
				inf = "{},{},{},{},{},{},{},{},{},{},{}:\n".format (k[0], k[1], base[k], depth, Mat, Mis, num_ins, Del, q_lst, k[2], bases_counts)
				outh.write (inf)
			except:
				sys.stderr.write ("file {} {} does not work\n".format (tsv,k))

def df_is_not_empty(df):
	'''
	input df is a df filtred on reference id 
	if is is empty: next (df.iterrows()) does not work
	otherwise it returns a row of df 
	'''
	try:
		next (df.iterrows())
		return True
	except:
		return False
	
def df_proc (df, outfn):
	'''
	input is a dataframe for either forward or reverse strand
	'''
	if not df_is_not_empty (df):
		print ("empty dataframe for {}".format(outfn), file=sys.stderr)
		return None
	outfh = open (outfn, 'w')
	header = "#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del"
	print (header, file=outfh)
	gb = df.groupby(['#Ref','pos','base','strand']).agg({
               'cov':['sum'],
               'mis':['sum'],
               'ins':['sum'],
               'del':['sum'],
               'qual':['sum']})
	gb = gb.reset_index()
	for row in gb.itertuples():
		coor = ",".join (map (str,row[1:5]))
		cov, mis, ins, _del, qual = row[5:]
		mis = '%0.5f' % (mis/cov)
		ins = '%0.5f' % (ins/cov)
		_de = '%0.5f' % (_del/cov)
		q = np.array ([x for x in qual.split(':') if x ]).astype(int)  
		qmn,qme,qst = '%0.5f' % np.mean(q), '%0.5f' % np.median(q), '%0.5f' % np.std(q)	
		outfh.write ("{},{},{},{},{},{},{},{}\n".format(coor,cov,qmn,qme,qst,mis,ins,_de))

def tsv_generator (reference, bam_file, sam2tsv, mapping_type):
	cmds = java_bam_to_tsv (bam_file, reference, sam2tsv, mapping_type)
	tsv_gen = []
	if mapping_type == 't':
		cmd = subprocess.Popen ((cmds[0]), stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True )
		returncode = cmd.returncode
		if returncode:
			print (res[1], file=sys.stderr)
			exit()
		tsv_gen.append (stdin_stdout_gen (cmd.stdout))
	elif mapping_type == 'g':
		run1 = subprocess.Popen ((cmds[0]), stdout=subprocess.PIPE, stderr = subprocess.PIPE,shell=True)
		run2 = subprocess.Popen ((cmds[1]), stdout=subprocess.PIPE, stderr = subprocess.PIPE,shell=True)
		if any ([run1.returncode, run2.returncode]):
			res1 = run1.communicate()
			res2 = run2.communicate()
			print (res1[1], res2[1], file=sys.stderr)
			exit()
		tsv_gen.append (stdin_stdout_gen(run1.stdout))
		tsv_gen.append (stdin_stdout_gen(run2.stdout))
		#tsv_gen = itertools.chain (stdin_stdout_gen (run1.stdout), stdin_stdout_gen (run2.stdout))
	return tsv_gen
	
def _prepare_dir (dirname):			
	if  os.path.exists(dirname):
		shutil.rmtree (dirname)
		sys.stderr.write ("{} already exists, will overwrite it\n".format(dirname))
		os.mkdir (dirname)
	if not os.path.exists (dirname):
		os.mkdir (dirname)

#~~~~~~~~~~~~~~~~~~~~~~~ main () ~~~~~~~~~~~~~~~~~~~~~~~
def main ():
	parser = argparse.ArgumentParser()
	required_args = parser.add_argument_group ('Required Arguments')
	required_args.add_argument ('-R','--reference', help='''samtools faidx indexed reference file and with
                        sequence dictionary created using picard
                        CreateSequenceDictionary''')
	required_args.add_argument ('-b', '--bam', type=str, help='bam file; if given; no need to offer reads file; mapping will be skipped')
	required_args.add_argument ('-s', '--sam2tsv',type=str, default='',help='/path/to/sam2tsv.jar; needed unless a sam2tsv.jar produced file is already given')
	parser.add_argument ('-n', '--number_cpus', type=int, default=4, help='number of CPUs') 
	parser.add_argument ('-T', '--type', type=str, default="t", help="reference types, which is either g(enome) or t(ranscriptome);")
	args=parser.parse_args()

#~~~~~~~~~~~~~~~~~~~~~~~ prepare for analysis ~~~~~~~~~~~~~~ 
	tsv_gen = None  # generator 
	prefix = '' 
 
	if args.reference:
		if not file_exist (args.reference):
			sys.stderr.write (args.reference + ' does not exist\n')
			exit()
		ref_faidx = args.reference +'.fai'
		if not file_exist (ref_faidx):
			sys.stderr.write (ref_faidx + ' needs to be created with samtools faidx\n')
			exit()	
		dict_fn = args.reference + '.dict'
		if not file_exist (dict_fn):
			sys.stderr.write (dict_fn + ' needs to be created using picard.jar CreateSequenceDictionary\n')
			exit()	
	if args.bam:
		bam_file = args.bam 
		if not file_exist (bam_file):
			sys.stderr.write (bam_file+' does not exist; please double check!\n')
			exit()
		else:
			if not file_exist (args.sam2tsv):
				sys.stderr.write ("Please offer correctly path to sam2tsv.jar\n".format(args.sam2tsv))
				exit()
			if not os.path.exists (bam_file+'.bai'):
				sys.stderr.write ('bam file not indexed!\nstarting indexing it ...')
				os.system ('samtools index ' + bam_file + '.bai')
			if not args.reference :
				sys.stderr.write('requires reference file that was used for reads mapping\n')
	if args.sam2tsv:
		sam2tsv = args.sam2tsv
		if not file_exist (sam2tsv):
			print (sam2tsv,'does not exist, Please provide it properly.', file=sys.stderr)
			exit()


	prefix = args.bam.replace ('.bam','')
	tsviter_lst = tsv_generator (args.reference, args.bam, args.sam2tsv, args.type)

	strands = ["plus_strand", "minus_strand"]

	tmp_dir = prefix + '_TMP_'
	_prepare_dir (tmp_dir)

	def tsv_to_var (tsvit, tmp_dir, out_var_fn, number_threads):
		processes = [] 	
		manager = Manager()
		q = manager.Queue(number_threads)
		ps = Process (target = split_tsv_for_per_site_var_freq, args = (tsvit, tmp_dir, q, number_threads, 2000))
		processes.append (ps) 
		for _ in range(number_threads):
			ps = Process (target= tsv_to_freq_multiprocessing_with_manager, args = (q, tmp_dir))
			processes.append (ps)
		for ps in processes:
			ps.daemon = True
			ps.start()
		for ps in processes:
			ps.join()		
		df = dd.read_csv ("{}/small_*freq".format(tmp_dir))
		out = out_var_fn 
		df_proc (df, out_var_fn)
		if  os.path.exists(tmp_dir):
			pool = mp.Pool(number_threads)
			tmp_files = glob.glob("{}/small*".format(tmp_dir))
			pool.map(_rm,  tmp_files)
			shutil.rmtree(tmp_dir)
	
	for idx, tsvit in enumerate (tsviter_lst):
		tmp_dir = prefix + '_TMP_'
		_prepare_dir (tmp_dir)
		out_var_fn = "{}.{}.per.site.csv".format (prefix, strands[idx])
		tsv_to_var (tsvit,  tmp_dir, out_var_fn, args.number_cpus)
	
	
if __name__ == "__main__":
	main()
