#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np 
import re

import sys

usage = "_USAGE_:	python {} {} {} {}".format ("Epinano_delta_sumErr.py", "WT_sumErr", "KO_sum_Err", "output")
usage += """
		applicable for results from Epinano_sum_err.py either single sites or kmer wsie...

"""

if len (sys.argv) < 4:
	print (usage, file=sys.stderr)
	exit (0)

def openfile (fn):
	if fn.endswith (".gz"):
		df = pd.read_csv (fn,compression='gzip')
	else:
		df = pd.read_csv (fn)
	return df 

# 5mer sumerr
#Ref,pos,base,strand,cov,sum_err1,sum_err2,sum_err3,sum_err4,sum_err5

ko=sys.argv[2]
wt=sys.argv[1]
out = sys.argv[3]


ko_df, wt_df = openfile (ko), openfile(wt)   			#pd.DataFrame ([]), pd.DataFrame([])
ko_df = ko_df.drop(['cov'], axis=1)
wt_df = wt_df.drop(['cov'], axis=1)

kdf = ko_df.set_index (["#Ref", 'pos', 'base','strand'])
wdf = wt_df.set_index (["#Ref", 'pos', 'base','strand'])

sub = wdf.sub (kdf)
sub = sub.reset_index()
sub = sub.dropna (how='any',axis='rows')
p = re.compile (r'sum_err(\d*)')
for col in sub.columns:
	m = p.match (col)
	if m:
		number = m.groups()[0]
		sub = sub.rename (columns={col:"delta_sum_err{}".format(number)})

		
#sub = sub.rename (columns={"sum_err":"delta_sum_err"})
#sub = sub.sort_values (by = 'delta_sum_err' , ascending=False)
sub.to_csv (out, index=False)
