#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from epinano_modules import openfile 
import matplotlib.pyplot as plt 
plt.style.use('ggplot')
import re, argparse, time 
from itertools import zip_longest, count
import scipy.stats as stats 
#from scipy.stats import ks_2samp
#from scipy.stats import mannwhitneyu
#from scipy.stats import kruskal
import numpy as np
#from scipy.stats import rankdata
import re 

#sys.setrecursionlimit(1000)
def statistical_methods_available():
	stats_methods  = """
	The Mann-Whitney U test for comparing independent data samples: 
		the nonparametric version of the Student t-test.
	The Wilcoxon signed-rank test for comparing paired data samples: 
		the nonparametric version of the paired Student t-test.
	The Kruskal-Wallis H and Friedman tests for comparing more than two data samples: 
		the nonparametric version of the ANOVA and repeated measures ANOVA tests.
	"""
	print (stats_methods, file = sys.stderr)


def clean_squiggle (s):
	ary = [i for i in s.split(',') if i and re.match(r'\d+', i)]
	return np.array(ary).astype(float)

def standard_two_samples_ttest (list1, list2):
	stat = stats.ttest_ind(list1, list2)
	return stat.pvalue

def mwu (s1, s2):
	stat = stats.mannwhitneyu (s1, s2)
	p1, p2 =  stat.statistic,  stat.pvalue
	return p2
	#return '{}\t{}'.format (p1, p2)

def ks2samp (s1, s2):
	stat = stats.ks_2samp (s1, s2)
	return  stat.pvalue
	#return '{}\t{}'.format (p1, p2)

def fdr(p_vals):
	'''
	p_vals: list of p-values 
	'''
	p_filt = np.array ([p for p in p_vals if p is not 'NaN'])
	ranked_p_values = stats.rankdata(p_filt)
	fdr = p_filt * len(p_vals) / ranked_p_values
	fdr[fdr > 1] = 1	
	fdr = (p for p in fdr)
	fdr_with_na = (next(fdr) if p is not 'NaN'  else 'NaN' for p in p_vals)	
	return fdr_with_na 

def homogenity_check (list_of_samples):
	'''
	[[], [], []] contains at least two list of intensity / dwell time values 
	'''
	pvalue = stats.levene (*list_of_samples)
	return pvalue >0.05 
	
def normality_check (list_of_values):
	'''
	[] contain data points of a sample subject to be tested 
	'''
	pvalue = stats.shapiro (list_of_values)
	return pvalue >0.05 
	
def median_test (list_of_samples):
	'''
	[[],[],..]
	'''
	stat, p, med, tbl = stats.median_test(*list_of_samples)
	return p
	
def krustal_test (list_of_samples):
	'''
	[[],[],...]
	'''
	stat = stats.kruskal (*list_of_samples)
	return stat.pvalue
	
def remove_outlier(list_of_values):
	z = np.abs (stats.zscore (list_of_values))
	outlier_idx = np.where (z>3)[0]

def main():
	parser = argparse.ArgumentParser("Stattistically compare >=2 samples in terms of current intensity valeus/dwell timedistributions;")
	parser.add_argument('--group1_files', type=str, required =True, 
									help='a series of samples of the same type, e.g. modified or unmodified; multiuple samples have to be seperated with comma')
	parser.add_argument ('--group2_files',type=str, required=True, help="In contrast to sample_group1, it accepts the oppositie type of samples ")
	#parser.add_argument ('--label1', type=str, required=True, help='laeble for group1 files')
	#parser.add_argument ('--label2', type=str, required=True, help='laeble for group2 files')
	parser.add_argument ('--output', type=str,required=True,help='output file')
	args = parser.parse_args()
	
	raw_values = args.output+'.raw'
	
	ks_c, ks_t, mwu_c, mwu_t, median_c, krustal_c, median_t, krustal_t = [], [] ,[],[],[],[],[],[]
	
	ref_sites = []
	
	group1, group2 = args.group1_files, args.group2_files
	group1 = re.sub (r'^,+','',group1)
	group1 = re.sub (r',+$','',group1)
	group2 = re.sub (r'^,+','',group2)
	group2 = re.sub (r',+$','',group2)
	
	group1, group2 = group1.split(','), group2.split(',')
	
	stats_methods = []
	if len(group1) == len(group2) == 1:
		# you can perform non-parameteric test of medians with ks2samp and mwu 
		stats_method = [ks2samp, mwu, krustal_test, median_test]
	else:
		stats_method = [krustal_test, median_test]
	
	raw = open (raw_values,'w')
	print ('Ref\tKs2samp_Intensity\tMWU_Intensity_Pvalue\tMedian_test_Intensity\tKrsutal_test_Intensity\t'
			'Ks2samp_Duration\tMWU_Duration\tMedian_test_Duration\tKrusta_test_Duration', file= raw)
	
	fhs = [openfile (f) for f in group1 + group2]
	for lines in zip_longest (*fhs, fillvalue=''):
		if lines[0].startswith ('index_'):
			continue
		group1_effective_lines =np.sum ([len(line.rstrip().split()) >1  for line in lines[:len(group1)]])
		group2_effective_lines =np.sum ([len(line.rstrip().split()) >1  for line in lines[len(group1):]])
		idx = lines[0].rstrip().split()[0]
		#not_enough_cols =  any ([len(line.rstrip().split()) ==1 for line in lines])
		if group1_effective_lines ==0: 
			print (f'warning{idx}does not have info in all {group1}', file=sys.stderr)
			continue
		elif group2_effective_lines == 0:
			print (f'warning{idx}does not have info in all {group2}', file=sys.stderr)
			continue
		cp1, cp2, cp3,cp4, dp1, dp2,dp3,dp4 = '','','','','','','',''
		if group1_effective_lines == 1 and group2_effective_lines == 1:
			line1 = lines[0].rstrip().split()
			c1, d1 = line1[-2], line1[-1]
			line2 = lines[1].rstrip().split()
			c2,d2 = line2[-2], line2[-1]
			d1, d2, c1, c2 = clean_squiggle(d1), clean_squiggle(d2), clean_squiggle(c1), clean_squiggle(c2)
			try:
				cp1 = ks2samp(c1,c2)
			except:
				print (line1[0], line2[1], 'current ks2sam failed', file=sys.stderr)
			try:
				dp1 = ks2samp (d1, d2)
			except:
				dp1 = 'NaN'
				print (line1[0], line2[1], 'duration ks2sam failed', file=sys.stderr)
			try:
				cp2 = mwu (c1,c2)
			except:
				cp2 = 'NaN'
				print (line1[0], line2[1], 'current mannwh test failed', file=sys.stderr)
			try:
				dp2 = mwu(d1, d2)
			except:
				dp2 = 'NaN'
				print (line1[0], line2[1], 'duration mannwh test failed', file=sys.stderr)
		else:
			cp1, dp1, cp2, dp2 = "NaN" , 'NaN' , 'NaN',  'NaN'
		
		try:
			cp3 = median_test([c1,c2])
		except:
			cp3 = "NaN"
			print (line1[0], 'current median test failed', file=sys.stderr)
		try:
			dp3 = median_test([d1, d2])
		except: 
			dp3 = 'NaN'
			print (line1[0], 'duration median test failed', file=sys.stderr)
		try:
			cp4= krustal_test ([c1,c2])
		except:
			cp4 = 'NaN'
			print (line1[0], 'current krustal_test failed', file=sys.stderr)
		try:
			dp4  = krustal_test([d1, d2])
		except:
			dp4 = 'NaN'
			print (line1[0], 'duration krustal_test failed', file=sys.stderr)
		del line1
		del line2 
		ks_c.append (cp1); ks_t.append(dp1)
		mwu_c.append(cp2); mwu_t.append(dp2)
		median_c.append (cp3); median_t.append(dp3)
		krustal_c.append(cp4); krustal_t.append(dp4)
		print (idx, cp1, cp2, cp3, cp4, dp1, dp2, dp3, dp4, sep = "\t", file=raw)		
		ref_sites.append(idx)
		
	raw.close()
		
	ks_c_fdr = fdr (np.array (ks_c).astype(float))
	ks_d_fdr = fdr (np.array (ks_t).astype(float))
	mwu_c_fdr = fdr (np.array(mwu_c).astype(float))
	mwu_d_fdr = fdr (np.array (mwu_t).astype(float))
	median_c_fdr = fdr (np.array(median_c).astype(float))
	median_t_fdr = fdr (np.array(median_t).astype(float))
	krustal_c_fdr = fdr (np.array(krustal_c).astype(float))
	krustal_t_fdr = fdr (np.array(krustal_t).astype(float))
	
	out = open (args.output,'w') 
	print ('Ref\tKs2samp_Intensity\tMWU_Intensity_Pvalue\tMedian_test_Intensity\tKrsutal_test_Intensity\t'
			'Ks2samp_Duration\tMWU_Duration\tMedian_test_Duration\tKrusta_test_Duration', file= out)
	for i, j in enumerate (ref_sites):
		print (j , next (ks_c_fdr),  next (mwu_c_fdr), next(median_c_fdr), next(krustal_c_fdr),  
				next(ks_d_fdr), next(mwu_d_fdr), next(median_t_fdr), next(krustal_t_fdr), sep="\t", file=out)
					
if __name__ == '__main__':
	main()
	










