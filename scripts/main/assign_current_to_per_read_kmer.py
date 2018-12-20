import sys
mem = {}

read_pos_intensity = sys.argv[1]
with open (sys.argv[1],'r') as fh:
	for  l in fh:
		if l.startswith ('#'):
			continue
		ary = l.strip().split(',')
		rd, pos, current = ary
		k = rd+':'+pos 
		mem[k] = current

header = '#ReadKmer,RefKmer,Ref,RefPos,Read,ReadPos,Q,M,I,D,CurrentIntensity'
print header
single_read_var_slided = sys.argv[2]
with open (sys.argv[2],'r') as fh:
	for l in fh:
		if l.startswith ('#'):
			continue
		ary = l.strip().split(',')
		rd, pos = ary[4],ary[5]
		pos = ":".join (pos.split('|'))
		k = rd+':'+pos
		if  k in mem:
			print l.strip()+','+mem[k] #ary[0]+'\t'+mem[ary[0]]
		else :
			print l.strip()+','+'NA'
