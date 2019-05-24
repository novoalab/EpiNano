import sys
import gzip
mem = {}

usage = " ". join(["python", sys.argv[0] , "read-kmer-current-intensity", "per read variation", ">output"])

def fopen (f):
    if f.endswith ('.gz'):
	fh = gzip.open (f,'r')
    else:
	fh = open (f,'r')
    return fh

if len (sys.argv) < 2:
    print usage
    exit()

read_pos_intensity = sys.argv[1]
fh = fopen (read_pos_intensity)

#with open (sys.argv[1],'r') as fh:
if (True):
	for  l in fh:
		if l.startswith ('#'):
			continue
		ary = l.strip().split(',')
		rd, pos, current = ary[:3]
		k = rd+':'+pos 
		mem[k] = current

fh.close()
header = '#Read,Read_Window,ReadKmer,Ref,RefKmer,Ref_Window,Q,M,I,D,CurrentIntensity'
print header
single_read_var_slided = sys.argv[2]
fh = fopen (single_read_var_slided)
#with open (sys.argv[2],'r') as fh:
if (True):
	for l in fh:
		if l.startswith ('#'):
			continue
		ary = l.strip().split(',')
		rd, pos = ary[0],ary[1]
		#pos = ":".join (pos.split('|'))
		k = rd+':'+pos
		if  k in mem:
			print l.strip()+','+mem[k] #ary[0]+'\t'+mem[ary[0]]
		else :
			print l.strip()+','+'NA'
fh.close()
