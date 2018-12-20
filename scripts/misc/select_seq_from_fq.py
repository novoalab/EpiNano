import sys

mem = {}

with open (sys.argv[1],'r') as fh:
	for l in fh:
		ary = l.strip().split()
		mem[ary[0]] = True

fh = open (sys.argv[2],'r')
for l in fh:
	h = l.strip()			
	rd = h.split()[0].replace('@','')
	s = fh.next().strip()
	n = fh.next().strip()
	q = fh.next().strip()
	if rd in mem:
		print h
		print s
		print n 
		print q
	
