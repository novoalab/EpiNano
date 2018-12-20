import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re,sys
import seaborn as sns
import math

#header
#ref	pos	refBase	sample1_base	sample1_base_qualities	sample2_bases	sample2_base_qualities	U	p_value	10:mean1	11:mean2	12:median1	13:median2	14:num_insertions1	15:num_deletions2	16:num_mism1	17:num_insertions2	18:num_deletions2	19:num_misms2

'''
mismatches
insertions
deletions
mean
median
will be plotted
'''
out = sys.argv[2]
outh = open (out,"w")

print >>outh , "samples\tq_mean\tq_median\tmismatches\tinsertions\tdeletions"


def log (v):
    if v == 0:
        return '0'
    else:
        return str (math.log(v))

with open (sys.argv[1],'r') as fh:
    for line in fh:
        if not re.search (r'^chr\d+',line):
            line.rstrip()
            eles = line.split ()
            depth1 = len (eles[3].split(',') )
            depth2 = len (eles[6].split(',') )
            eles[15] = log(float (eles[15]))
            eles[13] = log(float (eles[13]))
            eles[14] = log(float (eles[14]))
            eles[18] = log(float (eles[18]))
            eles[16] = log(float (eles[16]))
            eles[17] = log(float (eles[17]))

            print >>outh, "\t".join(["wt",eles[9],eles[11],eles[15],eles[13],eles[14]])
            print >>outh, "\t".join(["ko",eles[10],eles[12],eles[18],eles[16],eles[17]])

fh.close()
outh.close ()
'''
violin plot
df = pd.read_table (sys.argv[2])
fontsize = 10
fig, axes = plt.subplots()
header = list(df)
del header[0]

fontsize = 10
fig, axes = plt.subplots()
sns.violinplot(header[0],'deletions',data=df,ax=axes)
axes.set_xlabel(header[0])
axes.set_ylabel('deletions')
plt.savefig ("5.pdf",format="pdf")
'''

df = pd.read_table (sys.argv[2])
df1 = df.select_dtypes([np.int, np.float])

h = list (df)
tmp = sys.argv[2]
types = ['mean','median','mismatches','insertions','deletions']
for i, col in enumerate (df1):
     fontsize = 10
     fig, axes = plt.subplots()
     sns.violinplot(h[0],col,data=df,ax=axes)
     axes.set_xlabel(h[0])
     axes.yaxis.grid(True)
     axes.set_ylabel(col)
     out = tmp+'.'+types[i]+".png"
     plt.savefig (out,format="png")

