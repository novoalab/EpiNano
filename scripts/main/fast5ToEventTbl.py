import sys
try:
	import h5py
except:
	print "unable to import h5py, please install it"
	exit(1)
if len (sys.argv) < 2:
	print 'Usage: python2.7',sys.argv[0],'read.fast5', '> read.event.tbl'
	exit(1)

def get_readid (h5):
    k = h5['Raw/Reads/'].keys()[0]
    rid = h5['Raw/Reads/'+str(k)].attrs['read_id']
    return rid


def get_base_call_path (h5):
    base_call = ''
    for i in h5['Analyses/'].keys():
        if i.startswith('Basecall_'):
            base_call = i
        return base_call


def main ():
    inp = sys.argv[1]
    f5 = h5py.File(inp,'r')
    read_id = get_readid(f5)
    base_call = get_base_call_path(f5)
    if f5['Analyses/'+base_call+'/BaseCalled_template/Fastq'].value.startswith('@'):
        header = f5['Analyses/'+base_call+'/BaseCalled_template/Events'].value.dtype.names
        header =  'read_id\t'+'\t'.join (list (header))
        print header
        for i in f5['Analyses/'+base_call+'/BaseCalled_template/Events'].value:
            i = str(i).replace("(","")
            i = str(i).replace(",","\t")
            i = i[:-1]
            print read_id+'\t'+i

if __name__  == '__main__':
    main()



