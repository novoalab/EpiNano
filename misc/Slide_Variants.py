#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import sys, os 
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from epinano_modules import slide_per_site_var 

usage = '''
python Slide_Variants.py per_site_var kmer_length

'''

if (len(sys.argv) < 3):
	print (usage)
	print ("please provide 1) variants table from Epinano_Variants and 2) windown size(integer)\n", file=sys.stderr)
	sys.exit(0)

per_site_var = sys.argv[1]
kmer_len = int(sys.argv[2])

slide_per_site_var(per_site_var, kmer_len)


