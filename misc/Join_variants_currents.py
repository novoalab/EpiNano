
#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import pandas as pd 
import sys, argparse, fileinput, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from epinano_modules import openfile 
'''
combine variants feautures and squiggle features
q mis ins del 
median current, median duration 
StdCurrent, StdDuration 
'''

parser = argparse.ArgumentParser('''combine variants feautures and squiggle features
		q, mis, ins, del, median current, median duration 
		StdCurrent, StdDuration ''')
parser.add_argument('--variants',required=True, type=str, help='variants table')
parser.add_argument('--intensity',nargs='*', metavar='FILE', 
						help='file or piped stream of collapsed nanopolish eventalign results meidan and std of current intensities and durations')
parser.add_argument ('--outfile',default='combined feature table name', type=str, required = True)
args = parser.parse_args() 

if isinstance (args.intensity, str):
	if not os.path.exists (args.intensity):
		print (f"{args.infile} does not exist", file=sys.stderr)
		exit()

intensityfh = fileinput.input(files=args.intensity if args.intensity else ('-', ))
#if not args.intensity:
#	infilefh = os.fdopen(os.dup(sys.stdin.fileno()))

var = dict() # store variants 
cur = dict() # store currnet intensities; current intensity and duration table include all sites 
mem = dict()
outfh = open (args.outfile,'w')
with openfile (args.variants) as variant:
	for line in variant:
		if line.startswith ('#'):
			line = line.rstrip()
			print (line,'I1,I2,I3,I4,I5,D1,D2,D3,D4,D5', sep =',',file=outfh)
			continue 
		ary = line.rstrip().split(',')
		win,ref,strand = ary[1:4] # 1-based position
		var[(win, ref, strand)] = line.rstrip()
for line in intensityfh:
	if line.startswith ('#'):
		continue 
	ary = line.rstrip().split(',')
	i1,i2,i3,i4,i5,d1,d2,d3,d4,d5 = ary[4:14]
	try:
		win,ref,strand = ary[1:4] 
		begin, end = win.split('-')
		begin = str (int (begin) + 1)
		end = str (int (end) + 1)   ### 0-based --> 1-based position
		win = "{}-{}".format(begin, end)
	except:
		print (line.rstrip())
		
	k = (win,ref,strand)
	if k in var:
		print (var[k], ",".join (ary[4:]), sep=',', file = outfh)
