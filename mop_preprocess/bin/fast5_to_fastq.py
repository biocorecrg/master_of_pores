#!/usr/bin/env python
import sys, warnings,os
import numpy as np
import re
import gzip
import h5py, ont_fast5_api
from ont_fast5_api.fast5_interface import get_fast5_file

usage = '''
_EUSAGE_: python fast5_to_fastq.py fast5file
_INPUT_: both single- and multi- fast5 files can be processed
_BUGS_: report to elzedliu@gmail.com
'''

if len(sys.argv) < 2:
    print (usage); exit(0)

f5file = sys.argv[1]
hdf5 = h5py.File(f5file,'r')
f5 = get_fast5_file (f5file,mode='r')
reads_in_f5 = f5.get_read_ids()
number_of_reads =  len( reads_in_f5)


fq_out = open ('.'.join(f5file.split('.')[:-1])+'.fastq','w')

def fastq_from_fast5 (h5f):
    '''
    for single fast5 file, it is hdf5 file open with h5py
    for multi- fast5 file, it is hdf5['read_id'] for single read
    '''
    k = 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
    if k in h5f:
        fastq = h5f[k].value.decode()
        return fastq
    else:
        return None

if number_of_reads == 1:
    fq = fastq_from_fast5 (hdf5)
    if fq :
        fq_out.write (fq)
    else:
        sys.stderr.write("no sequence basecalled for {} \n".format(f5file))
elif number_of_reads > 1:
    for rd in hdf5.keys():
        fq = fastq_from_fast5 (hdf5[rd])
        if fq:
            fq = fq.replace('_Basecall_1D_template','')
            fq_out.write(fq)
        else:
            sys.stderr.write("no basecall for read {} \n" .format(rd))
fq_out.close()
