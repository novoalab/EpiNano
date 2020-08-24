#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import sys, os, inspect
from inspect import getsourcefile  
from os.path import abspath  
from collections import defaultdict 
import scipy.stats as stats
import numpy as np

#currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(currentdir)
#src_dir = abspath(getsourcefile(lambda:0))
parent_dir = os.path.dirname (os.path.dirname (abspath(sys.argv[0])))
sys.path.insert(0, parent_dir)
from epinano_modules import openfile

mod = sys.argv[1] 
unm = sys.argv[2]
cov = 5

if len (sys.argv) == 4:
	try:
		cov = int(sys.argv[3])
	except:
		raise
		sys.exit(1)

mem = dict() #defautldict(list)
##Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del
with openfile (mod) as fh:
	for l in fh:
		if l.startswith ('#'):
			continue
		ary = l.rstrip().split(',')
		if int(ary[4]) < cov:
			continue
		idx = tuple (ary[:4])
		ary[10] = float (ary[10]) / 2
		err =  np.sum ([round (float(i) * int(ary[4])) for i in ary[8:10]])
		mem[idx] = [int(ary[4])] +[err]  

print ("#Ref,pos,base,strand,FET_Pvalue")
with openfile (unm) as fh:
	for l in fh:
		if l.startswith ('#'):
			continue 
		ary = l.rstrip().split(',')
		cov2 = int (ary[4])
		if cov2 < cov:
			continue
		idx = tuple (ary[:4])
		ary[10] = float(ary[10]) / 2
		if idx in mem:
			cov1, err1 = mem[idx]
			correct1 = cov1 - err1 #mem[idx][1], mem[idx][0] - mem[idx];  mem[idx][0], np.array(mem[idx][1:]) #mem[idx][0], mem[idx][1], mem[idx][2], mem[idx][3]
			err2 = np.sum ([round (float(i) * cov2) for i in ary[8:]])
			correct2 = cov2 - err2
			#print (cov1, cov2, err1, err2, correct1, correct2)
			try:	
				oddsratio, pv = stats.fisher_exact ([[err1, correct1], [err2, correct2]])
				del mem[idx]
				print (",".join (list(idx)), '%.3f' % pv, sep=",")
			except:
				print ('analysis of', idx, "failed", file=sys.stderr)
