#!/usr/bin/env python 
import sys, os,gzip,bz2,argparse 
from collections import OrderedDict

def openfile(f):
    if f.endswith ('.gz'):
        fh = gzip.open (f,'rt')
    elif f.endswith ('bz') or f.endswith ('bz2'):
        fh = bz2.open(f,'rt')
    else:
        fh = open(f,'rt')
    return fh 

def print_last_consecutive_lines (lines, outfh):
	contents = OrderedDict()
	for line in lines:
		ary = line.strip().split(',')
		ary = line.strip().split()
		bry = ary[0].split(';')
		ref,pos,base,strand = bry[0], bry[1], bry[2][0], bry[3]
		if len (ary)>1:
			contents[(ref,pos,strand)] = ",".join ([ref, pos, base, strand, ary[3], ary[6]])
		else:
			contents[(ref,pos,strand)] = ",".join ([ref, pos, base, strand, 'NaN', 'NaN'])
	
	win = len (lines)
	ary = lines [win//2].rstrip().split()
	bry = ary[0].split(';')
	middle = [] 
	if len (ary) >1:
		middle = [bry[0], bry[1], bry[2][0],bry[3], ary[3], ary[6]]
	else:
		middle = [bry[0], bry[1], bry[2][0],bry[3],"NaN","NaN"]
	window = str(int(middle[1]) - win//2)+'-'+str(int(middle[1]) + win//2 )
	kmer = ''
	consecutive_lines = []
	ref,pos,base,strand = middle[:4]
	for i in reversed (list (range(1, win//2+1))):
		k = (ref,str(int(pos) -i), strand)
		relative_pos = '-'+str(i)
		if k in contents:
			kmer = kmer + contents[k].split(',')[2]
			consecutive_lines.append (window + ',' + relative_pos +','+contents[k])
		else:
			kmer = kmer + 'N'
			consecutive_lines.append (window +','+relative_pos+','+",".join (list(k)) + 'NaN,NaN')
	consecutive_lines.append (window +',' + '+0' +','+",".join (middle))
	kmer = kmer + middle[2]
	for i in range(1,win//2+1):
		k = (ref,str(int(pos) + i),strand)
		relative_pos = '+'+str(i)
		if k in contents:
			kmer = kmer + contents[k].split(',')[2]
			consecutive_lines.append (window + ',' + relative_pos+ ',' + contents[k])
		else:
			kmer = kmer + 'N'
			consecutive_lines.append (window +',' + relative_pos+','+",".join (list(k)) +',' + 'NaN,NaN')
	for l in consecutive_lines:
		print (kmer+','+l, file=outfh)
		

def slide_intensity (intensity_table_file ,win=5):
	'''
	index_column    model_kmer_counts        mean_intensity median_intensity        std_intensity    mean_duration  median_duration std_duration     intensity_samples      duration_samples
	'''
	prefix =  intensity_table_file   #+'.' + str(win)+'mers.csv'
	out_tmp = prefix +'.{}mer.tmp'.format(win)
	if os.path.exists (out_tmp):
		os.remove (out_tmp)
	outfh = open (out_tmp,'w')
	
	fh = openfile (intensity_table_file)
	
	fh.seek(0,os.SEEK_END)
	eof = fh.seek (fh.tell()-1, os.SEEK_SET)
	#print (eof, 'is end of file'); exit()
	fh = openfile (intensity_table_file)
	head = fh.readline () 
	
	lines = [] 
	for _ in range (win):
		l = fh.readline().rstrip()
		if l:
			lines.append (l)
	if len  (lines) < win:
		print ('not enough sites to be slided',file=sys.stderr)
	
	contents = OrderedDict()
	
	for line in lines:
		ary = line.rstrip().split()
		bry = ary[0].split(';')
		ref,pos,base,strand = bry[0], bry[1], bry[2][0], bry[3]
		if len (ary)>1:
			contents[(ref,pos,strand)] = [ref, pos, base, strand, ary[3], ary[6]]
		else:
			contents[(ref,pos,strand)] = [ref, pos, base, strand, 'NaN', 'NaN']
	 
	while (fh.tell() <= eof):
		kmer = ''
		ary = lines [win//2].split()
		bry = ary[0].split(';')
		middle = []
		if len (ary)>1:
			middle = [bry[0], bry[1], bry[2][0],bry[3], ary[3], ary[6]]
		else:
			middle = [bry[0], bry[1], bry[2][0],bry[3], "NaN", "NaN"]
		window = str(int(middle[1]) - win//2)+'-'+str(int(middle[1]) + win//2)
		consecutive_lines = []
		ref,pos,base,strand = middle[:4]
		k_to_del = (ref,str(int(pos)-win),strand)
		for i in reversed (list (range(1, win//2+1))):
			k = (ref,str(int(pos) -i),strand)
			relative_pos = '-'+str(i)
			if k in contents:
				kmer = kmer + contents[k][2]
				consecutive_lines.append (window+','+relative_pos +','+",".join (contents[k]))
			else:
				kmer = kmer +'N'
				consecutive_lines.append (window+','+relative_pos+','+ "," . join ([ref, str(int(pos) -i), 'N', strand,  'NaN,NaN']))
		consecutive_lines.append (window+',+0' +','+",".join (middle))
		kmer = kmer +middle[2]
		for i in range(1,win//2+1):
			k = (ref,str(int(pos) + i),strand)
			relative_pos = '+'+str(i)
			if k in contents:
				kmer = kmer +contents[k][2]
				consecutive_lines.append (window+','+relative_pos+','+",".join (contents[k]))
			else:
				kmer = kmer+'N'
				consecutive_lines.append (window+','+relative_pos+','+ "," . join ([ref, str(int(pos) +i), 'N', strand, "NaN,NaN"])) 
		for l in consecutive_lines:
			print (kmer+','+l,file = outfh)
		del consecutive_lines
		if k_to_del in contents:
			del contents[k_to_del]
		lines = lines[1:]
		new_line = fh.readline().rstrip() #
		lines.append (new_line)
		ary = new_line.split()
		bry = ary[0].split(';')

		ref,pos,base,strand = bry[0], bry[1], bry[2][0], bry[3]
		if len (ary)>1:
			contents[(ref,pos,strand)] = [ref, pos, base, strand, ary[3], ary[6]]
		else:
			contents[(ref,pos,strand)] = [ref, pos, base, strand, 'NaN', 'NaN']
	print_last_consecutive_lines (lines, outfh)
	outfh.close()
	
	out2 = prefix + '.{}mer.csv'.format(win)
	outh2 = open (out2,'w')
	intensity_in_head = ",".join (["E{}".format(i) for i in range(1,win+1)])
	duration_in_head = ",".join (["D{}".format(i) for i in range(1,win+1)])
	outh2.write ('#Kmer,Window,Ref,Strand,{},{}\n'.format(intensity_in_head, duration_in_head))
	# ~~~~~~~~~~~~~~~~~
	tmpfh = open (out_tmp,'r')
	firstline = tmpfh.readline().rstrip().split(',')
	current_win = "{},{},{},{}".format(firstline[0], firstline[1], firstline[3],  firstline[6])
	lines = []
	lines.append (firstline)
	for l in tmpfh:
		ary = l.rstrip().split(',')
		window = "{},{},{},{}".format(ary[0], ary[1], ary[3],ary[6])
		if window.startswith ('NNNNN'):
			continue 
		if window != current_win:
			intensity, duration  = [], []
			for ele in lines:
				intensity.append (ele[7])
				duration.append (ele[8])
			Is = ",".join (intensity)
			Ds = ",".join (duration)
			if ('NaN' not in Is and "NaN" not in Ds):
				print (current_win, Is, Ds, sep=",", file= outh2)
			current_win = window
			lines = []
		lines.append (ary)
	# last 5 lines
	Is,Ds  = [], []
	for ele in lines:
		intensity.append (ele[7])
		duration.append (ele[8])
	Is = ",".join (intensity)
	Ds = ",".join (duration)
	if (Is not in "NaN,NaN,NaN,NaN,NaN" and Ds not in "NaN,NaN,NaN,NaN,NaN" and ary[0] not in 'NNNNN'):
		print (window, Is, Ds, sep=",", file= outh2)
	
	tmpfh.close()
	outh2.close()
	os.remove (out_tmp)
	return (out2)

def main ():
	parser = argparse.ArgumentParser("create a table containing current intensity and duration information in 5mer windows")
	parser.add_argument ('--intensity_table', required=True, help="current intensity tab seperated file")
	parser.add_argument ('--window', type=int, default=5, help = "sliding window size, default is 5")
	args = parser.parse_args()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	inp = args.intensity_table
	win = args.window
	slide_intensity (inp, win)
if __name__ == '__main__':
	main()
	
	
