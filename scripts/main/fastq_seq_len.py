import sys
import gzip

if len (sys.argv) < 2:
	print 'Usage: python2.7', sys.argv[1],'in.fastq > in.seq.len'
	exit(1) 

fastq = sys.argv[1]
if fastq.endswith ('.gz'):
    fqh = gzip.open (fastq,'r')
else:
    fqh = open (fastq,'r')

for line in fqh:
    rd_id = line.strip().split()[0]
    rd_id = rd_id.replace('@','',1)
    seq = fqh.next().strip()
    _,_ = fqh.next(),fqh.next()
    print rd_id, len(seq)
