import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re,sys
import seaborn as sns
import math
'''
targetSite      q1      q2      mis1    mis2    del1    del2    ins1    ins2
'''
out = sys.argv[2]
outh = open (out,"w")
print >>outh , "samples\tquality\tmismatches\tdeletions\tinsertions"


def log (v):
    if v == 0:
        return '0'
    else:
        return str (math.log(v))

with open (sys.argv[1],'r') as fh:
    for line in fh:
        if re.search (r'^c',line):
            line.rstrip()
            eles = line.split ()
            depth1 = len (eles[3].split(',') )
            depth2 = len (eles[6].split(',') )
            q1 = float (eles[1])
            q2 = float (eles[2])
            mis1 = log(float (eles[3])+3)
            mis2 = log(float (eles[4])+3)
            del1 = log(float (eles[5])+3)
            del2 = log(float (eles[6])+3)
            ins1 = log(float (eles[7])+3)
            ins2 = log(float (eles[8])+3)

            print >>outh, "\t".join(map (str, ["mod",q1,mis1,del1,ins1]))
            print >>outh, "\t".join(map (str, ["unm",q2,mis2,del2,ins2]))

fh.close()
outh.close ()

df = pd.read_table (sys.argv[2])
df1 = df.select_dtypes([np.int, np.float])

h = list (df)
tmp = sys.argv[2]
types = ['quality','mismatches','deletions','insertions']
for i, col in enumerate (df1):
     fontsize = 10
     fig, axes = plt.subplots()
     sns.violinplot(h[0],col,data=df,ax=axes)
     axes.set_xlabel(h[0])
     axes.yaxis.grid(True)
     axes.set_ylabel(col)
     out = tmp+'.'+types[i]+".png"
     plt.savefig (out,format="png")

