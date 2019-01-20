import random
import sys
'''
def random_line(afile):
	line = next(afile)
    for num, aline in enumerate(afile, 2):
      if random.randrange(num): continue
      line = aline
    return line
'''
if (len (sys.argv) < 2):
	print "Usage: python", sys.argv[0] , "<input file> number_of_lines"

if len(sys.argv) == 3:
	num_of_lines =  int(sys.argv[2])
else:
	num_of_lines = 10000

def random_sampler(filename, k):
	sample = []
	with open(filename, 'rb') as f:
		head = f.readline() #skip header line
		print head,
		f.seek(0, 2)
		filesize = f.tell()
		random_set = sorted(random.sample(xrange(filesize), k))
		for i in xrange(k):
			f.seek(random_set[i])
			# Skip current line (because we might be in the middle of a line) 
			f.readline()
			# Append the next line to the sample set 
			sample.append(f.readline().rstrip())
	return sample

random_sample = random_sampler (sys.argv[1],num_of_lines)
for i in random_sample:
	print i
