#!/usr/bin/env python3
import warnings
import h5py
import sys
#with warnings.catch_warnings():#
#    warnings.filterwarnings("ignore",category=DeprecationWarning)
#   import h5py
warnings.filterwarnings("ignore",category=DeprecationWarning)
f5file = sys.argv[1]
hdf5 = h5py.File(f5file,'r')

if 'Analyses' in hdf5.keys (): #single read fast5 file
# <KeysViewHDF5 ['Analyses', 'PreviousReadInfo', 'Raw', 'UniqueGlobalKey']>
    try:
        k = list (hdf5['Raw/Reads/'].keys())[0]
        rid = hdf5['Raw/Reads/'+str(k)].attrs['read_id']
        if 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq' not in hdf5:
            sys.stderr.write(rid+' no basecall')
            exit(0)
        if 'Analyses/Basecall_1D_000/BaseCalled_template/Events' in hdf5:
            header = hdf5['Analyses/Basecall_1D_000/BaseCalled_template/Events'].value.dtype.names
            header =  '#read_id\t'+'\t'.join (list (header))
            events = hdf5['Analyses/Basecall_1D_000/BaseCalled_template/Events'].value
            print (header)
            for i in events:
                i = str(i).replace("(","")
                i = str(i).replace(",","\t")
                pat = """b'"""
                i = i.replace (pat,'')
                i = i.replace ("""'""",'')
                i = i.replace (""")""",'')
                print (str (rid.decode())+'\t'+i)
                print (str (rid.decode())+'\t'+i)
    except:
         sys.stderr.write(rid + ' event table extraction failed\n')
else:
    for rd in hdf5.keys():
        try:
            if 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq' not in hdf5[rd]:
                sys.stderr.write(rd+' no basecall')
                continue
            if 'Analyses/Basecall_1D_000/BaseCalled_template/Events' in hdf5[rd]:
                header = hdf5[rd]['Analyses/Basecall_1D_000/BaseCalled_template/Events'].value.dtype.names
                header =  '#read_id\t'+'\t'.join (list (header))
                events = hdf5[rd]['Analyses/Basecall_1D_000/BaseCalled_template/Events'].value
                print (header)
                for i in events:
                    i = str(i).replace("(","")
                    i = str(i).replace(",","\t")
                    pat = """b'"""
                    i = i.replace (pat,'')
                    i = i.replace ("""'""",'')
                    i = i.replace (""")""",'')
                    #i = i[:-1]
                    rd = rd.replace('read_','')
                    print (str (rd)+'\t'+i)
                    print (str (rd)+'\t'+i)
                sys.stderr.write(rd + ' event table extraction successful\n')
        except:
            raise
            sys.stderr.write(rd + ' event table extraction failed\n')
