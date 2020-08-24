#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import sys, glob, os 


def openfile(f):
    if f.endswith ('.gz'):
        fh = gzip.open (f,'rt')
    elif f.endswith ('bz') or f.endswith ('bz2'):
        fh = bz2.open(f,'rt')
    else:
        fh = open(f,'rt')
    return fh

def final_output (sumfiles, outfh):
	header = openfile (sumfiles[0]).readline().rstrip()
	print (header, file = outfh)
	for f in sumfiles:
		header = ''
		with openfile (f) as fh:
			for l in fh:
				if l.startswith('index_column'):
						continue 
				print (l.rstrip(),file=outfh)
		os.remove (f)
	outfh.close()	
	
def main ():
	input_dir  = sys.argv[1]
	outfh = open (f'{input_dir}/Intensity.collapsed.tsv','w')
	files = glob.glob (f'{input_dir}/*.sum')
	print ('combine', files)
	final_output (files, outfh)
if __name__ == '__main__':
	main()
