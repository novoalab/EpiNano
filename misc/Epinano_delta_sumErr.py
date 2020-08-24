#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np 

import sys

usage = "_USAGE_:	python {} {} {}".format ("Epinano_delta_sumErr.py", "WT_sumErr", "KO_sum_Err", "output")

if len (sys.argv) < 4:
	print (usage, file=sys.stderr)
	exit (0)

def openfile (fn):
	if fn.endswith (".gz"):
		df = pd.read_csv (fn,compression='gzip')
	else:
		df = pd.read_csv (fn)
	return df 

# when single sites
ko=sys.argv[1]
wt=sys.argv[2]
out = sys.argv[3]


ko_df, wt_df = openfile (ko), openfile(wt)   			#pd.DataFrame ([]), pd.DataFrame([])
ko_df = ko_df.drop(['cov'], axis=1)
wt_df = wt_df.drop(['cov'], axis=1)

kdf = ko_df.set_index (["#Ref", 'pos', 'base','strand'])
wdf = wt_df.set_index (["#Ref", 'pos', 'base','strand'])

sub = wdf.sub (kdf)
sub =sub.reset_index()
sub = sub.dropna (how='any',axis='rows')
sub = sub.rename (columns={"sum_err":"delta_sum_err"})
sub = sub.sort_values (by = 'delta_sum_err' , ascending=False)

sub.to_csv (out, index=False)
