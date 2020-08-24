#!/usr/bin/env python
from Bio import SeqIO
import sys, argparse,os,glob, fileinput, re 
from multiprocessing import Process, Queue, Pool
import dask
import dask.dataframe as dd
import numpy as np 
import pandas as pd 
from itertools import zip_longest 
import datetime,shutil 
from collections import Counter, defaultdict
# -*- coding: utf-8 -*-

# add py path
DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0,DIR)
from epinano_modules import openfile

#~~~~~~~~~~~~~~~~~~ Private functions ~~~~~~~~~~~~~~~~~
def _split_eventalign_tbl_on_read (fh , in_q, tmp_dir, strand, num_cpus=3, num_reads_in_each_chunk=4000):
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
		ary = fh.readline().rstrip().split()
		header = "\t".join ([ary[0], ary[1],ary[2],ary[4],ary[8],ary[9],ary[15]])
		print (header, file = outfh)
		for line in fh:
			ary = line.rstrip().split()
			if line and line.startswith ('contig'):
				continue
			rd = ary[3]
			ary[4] = strand 
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
		#print ('finish splitting input file', datetime.datetime.now(), file = sys.stderr)
		for i in range (num_cpus):
			in_q.put (None)

def load_single_file (file):
	'''
	load raw nanopolish event align output(s); usecols -->
	contig  position reference_kmer  strand model_kmer event_length samples
	'''
	columns ="contig  position reference_kmer  strand model_kmer event_length samples".split()
	if file.endswith ('.gz'):
		return (dd.read_csv (file, sep="\s+",usecols=columns, compression='gzip'))
	else:
		return (dd.read_csv (file, sep="\s+",usecols=columns))
		

def df_to_csv (df,outfile):
	try:
		df.to_csv (outfile,sep='\t', single_file=True, index=False, header=True)
	except:
		print (f'warning: df_to_csv --> {outfile} did not go well')
		print (f'warning: df_to_csv --> {outfile} did not go well')
		raise

def groupby (df):
	'''
	groubpy --> 'contig', 'position', 'reference_kmer',"strand" of nanopolish eventalign resutls 
	'''
	df.event_length = df.event_length.apply(str,meta='string')
	df.samples = df.samples.apply(str,meta='string')
	df.samples = df.samples +','	
	df.event_length = df.event_length+','
	df.model_kmer = df.model_kmer+','
	grp = df.groupby ([ 'contig', 'position', 'reference_kmer',"strand"]).agg({'samples':['sum'], 'event_length':['sum'], 'model_kmer':['sum']}).reset_index()
	grp.columns = grp.columns.droplevel(1)
	grp.position = grp.position.apply(str,meta=str)
	grp['index_column'] =  grp[['contig', 'position', 'reference_kmer', 'strand']].apply(lambda x: ';'.join(x),axis=1,meta=str)
	grp=grp.drop (['contig', 'position', 'reference_kmer', 'strand'],axis=1)
	return grp 
	
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

def add_ref_table_to_chunk_table (ref_file, ref_cols, ref_file_df, chunk_file, chunk_file_df):
	"""
	ref_file_df and chunk_file_df are indexed on the index_column
	ref_cosl is column with ref, pos, strand sorted 
	"""
	df1 = ref_file_df  #dd.read_csv(ref_file ,sep="\s+")
	df2 = chunk_file_df #dd.read_csv(chunk_file ,sep="\s+")
	seq_id = ref_file[:len(ref_file) - 10]  #.+.evn.tbl name offset is 10 -- suffix length 
	df2 = df2[df2.index_column.str.startswith(seq_id)]
	
	df1 = df1.set_index('index_column')
	df2 = df2.set_index('index_column')
	
	try:
		if df_is_not_empty(df2): # > 0: only process non-empty dataframe after filtering on ref id 
			#df1 = df1.applymap(str); #df2 = df2.applymap(str)
			sum = df1.add(df2)
			sum = sum.compute()
			sum = sum.reindex (ref_cols)
			sum = sum.reset_index()
			outfile = chunk_file+ '.' + ref_file 
			sum.to_csv (outfile, sep='\t', index=False, header=True)
		#else: #	print (f'{chunk_file} has no informtion relevant to {ref_file}', file=sys.stderr)
	except:
		print ('fail add ref table with chunk table',file = sys.stderr)
		raise 

def process_splitted_nanopolish_event_align_files (inq, initial_ref_tbls, num_cpus): #, matchq, mismatchq): #allq,matchq,mismatchq):
	f=''
	try: 
		for f in iter (inq.get,None):
			df = load_single_file(f)
			grp = groupby (df)
			grp_file = f + '.grp'
			try:
				for ref in initial_ref_tbls:
					ref_df = dd.read_csv(ref, sep="\s+")
					cols = ref_df['index_column'].compute()
					add_ref_table_to_chunk_table (ref, cols, ref_df, grp_file, grp)
			except:
				print (f'add {grp_file} and {initial_ref_tbls} failed', file=sys.stderr)
				raise
			os.remove (f)
	except:
		print (f'process single splitted nanopolish event align results file {f} failed', file = sys.stderr)
		raise 
				
def rmfile (fn):
	os.remove (fn)
	
def print_out (files_q,results_dir,sumfiles_q):
	for fn, files_lst in iter (files_q.get, None):
		sumfile = f'{results_dir}/{fn}.sum'
		filehandles = [openfile (f) for f in files_lst]
		outfh = open (sumfile,'w')
		print (
			('index_column\tmodel_kmer_counts\t'),
			('mean_intensity\tmedian_intensity\tstd_intensity\t'),
			('mean_duration\tmedian_duration\tstd_duration\t'),
			('intensity_samples\tduration_samples'), file = outfh)
		for lines in zip_longest (*filehandles, fillvalue=''):
			if lines[0].startswith('index_column'):
				continue 
			index_column, samples, event_length, model_kmer = '', [], [], []
			s, e = '',''
			info =defaultdict(list)
			info['model_kmer'] = []
			info['samples'] = []
			info['event_length'] = []
			for l in lines: 
				ary = l.rstrip().split('\t')
				index_column = ary[0]
				if len (ary)>1:
					s = s +',' + ary[1]
					e = e+','+ ary[2]
					info['model_kmer'] = info.get('model_kmer',[])+ary[3][4:].rstrip(',').split(',')
					info['samples'] = info.get('samples',[]) + ary[1][4:].rstrip(',').split(',')
					info['event_length'] = info.get('event_length',[])+ary[2][4:].rstrip(',').split(',')
		
			s = s.replace ('nan,','').rstrip(',') 
			e = e.replace ('nan,','').rstrip(',')
			samples = np.array (info['samples']).astype(float)
			event_length = np.array (info['event_length']).astype(float)
			if len (samples) == 0: # no line contains information 
				print (index_column,file = outfh) #'nan,','nan,','nan,','nan,','nan,','nan,','nan,','nan,','nan,',file=outfh)
			else:
				model_kmer = Counter (info['model_kmer'])	
				model_kmer  =",".join ([i+':'+str(j) for i, j in model_kmer.items()])
				print (index_column, model_kmer, np.mean(samples), np.median(samples), np.std(samples), 
								np.mean(event_length), np.median (event_length), np.std(event_length),
								s,e,sep="\t", file=outfh)
		outfh.close()
		for f in files_lst:
			os.remove (f)	
			#chk_file = os.path.basename(f).replace('.grp',' ').split()[0]
			#os.remove ("{}/{}".format (os.path.dirname(f), chk_file ))
			
		sumfiles_q.put(sumfile)
		
#def generate_final_output():
def initiate_ref_table (ref_fn, strand = '+'):
	initial_ref_tbls = []
	with openfile (ref_fn) as fh:
			for record in SeqIO.parse(fh,'fasta'):
				id,seq = record.id, record.seq
				initial_tbl = id+'.'+strand+'.evn.tbl'
				if not os.path.exists (initial_tbl):
					initial_table_fh = open (initial_tbl,'w')
					header = "index_column\tsamples\tevent_length\tmodel_kmer"
					print (header, file=initial_table_fh)
					for i in range (len(seq)):
						kmer=str (seq[i:i+5])
						idx = ";".join ([id,str(i),kmer,strand])
						print (idx,'nan,','nan,','nan,',sep="\t", file=initial_table_fh)
					initial_table_fh.close()	
				initial_ref_tbls.append (initial_tbl)
	return initial_ref_tbls
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def main ():
	parser = argparse.ArgumentParser('process nanopolish event align results for the purpose of analyzing them at per reference position level')
	parser.add_argument('--reference', required=True, type=str, help='reference fasta file used for nanopolish event align')
	#parser.add_argument('--infile', nargs="*", metavar='FILE', 
	#								help='input file (or linux pipe) of nanopolish eventalign results (with signal samples)')
	parser.add_argument ("--infile", type=str, help = "nanopolish eventalign results")
	parser.add_argument ('--outdir',  default='PerSiteIntensityAndDuration', type=str, required = True)
	parser.add_argument ('--chunking_size', default=4000, type = int, help = '''raw input will be splitted to multiple chunks and parallely processed;\n 
																													this argument controls the number of reads in one chunk ; default = 4000''')
	parser.add_argument ('--threads', type=int, default=3)
	#parser.add_argument ('--keep_tmp', action = 'store_true', help = 'keep intemediate data, only switch it on when you have enough disk space ')
	parser.add_argument ('--strand', type=str, default ='+', help = ('''the strand on which nanopolish eventalign results to keep;\n
														strand can be '+' or '-' or '+-' \n
				                                 		in the case of transcriptome mapping and eventalign, strand is '+';\n
														in the case of genome mppping, strand is '+-'; '''))
	args = parser.parse_args()
	
	print ('start analysis', datetime.datetime.now())
# ~~~~~~~~~~~~~~~~ process arguments ~~~~~~~~~~~~~~~~~~~~~~
	results_dir =  args.outdir+'_events.collapsed'
	tmp_dir = results_dir #args.outfile +'.events.collapsed'
	if not os.path.isdir (results_dir):
			os.mkdir (results_dir)
	else:
		print (f"{results_dir} exists, will over-write it", file=sys.stderr)
			
	if isinstance (args.infile, str):
		if not os.path.exists (args.infile):
			print (f"{args.infile} does not exist", file=sys.stderr )
			exit()
	
	infilefh = openfile (args.infile)
	"""
	infilefh = None 		
	if args.infile and args.infile[0].endswith('.gz'):
		infilefh = fileinput.input(files=args.infile )
	elif not args.infile[0].endswith ('.gz'):
		infilefh = fileinput.input(files=args.infile if args.infile else ('-', ))
	else:
		infilefh = os.fdopen(os.dup(sys.stdin.fileno()))
	"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   create single reference feature table 
	initial_ref_tbls = [] 
	for s in args.strand:	
		initial_ref_tbls = initial_ref_tbls + initiate_ref_table (args.reference, s)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	split nanopolish results on reads into small chunks 	
	procs = []
	inq = Queue(maxsize=args.threads)
	ps = Process (target = _split_eventalign_tbl_on_read, args = (infilefh, inq, tmp_dir, args.strand, args.threads, args.chunking_size))
	procs.append (ps)
	ps.start()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# apply groupby to each chunk, make it noredundant
# add each non-redundant small chunk all to single reference tables 

	print ('adding small chunk file and referenc file', datetime.datetime.now(),file=sys.stderr)
	for _ in range (args.threads):
		ps = Process (target = process_splitted_nanopolish_event_align_files, args = (inq, initial_ref_tbls, args.threads))
		ps.start()
		procs.append(ps)
	for ps in procs:
		ps.join()
	print ('finish adding reference and small files' , datetime.datetime.now(),file=sys.stderr)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	output single reference table 
	files_q = Queue()
	for fn in initial_ref_tbls:
		files = glob.glob (f"{tmp_dir}/*.{fn}")
		if files:
			files_q.put ((fn, files)) #(fn, zip_longest (*filehandles, fillvalue='')))
	for _ in range (args.threads):
		files_q.put(None)
	procs = []	
# output one reference one file 
	sumfiles_q = Queue()
	for _ in range (args.threads):
		ps = Process (target = print_out, args = (files_q,results_dir,sumfiles_q))
		ps.start()
	for ps in procs:
		ps.join()
		if ps.is_alive():
			print ('printing sum job is not finished')

# combine reference specific output to final output 
	#print ('finish outputting single reference associated event table',datetime.datetime.now(), file = sys.stderr)
	print ('finish warapping up all small files ', datetime.datetime.now(), file = sys.stderr)
	print ('if you want to CAT all results together, run concat_events.py to concatenate the results')
	

if __name__ == '__main__':
	main()





				
				


