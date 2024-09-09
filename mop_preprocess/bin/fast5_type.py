#!/usr/bin/env python
import os
import sys
import warnings

import h5py
import ont_fast5_api
from ont_fast5_api.fast5_interface import get_fast5_file

__author__ = 'Huanle.Liu@crg.eu'
__version__ = '0.1'
__email__ = 'same as author'

usage = '''
        python fast5_type.py fast5file
        return:
            0: single read fast5
            1: multi-reads fast5
'''

if len (sys.argv) < 2:
    print (usage)
    exit(0)

f5file = sys.argv[1]
hdf5 = h5py.File(f5file,'r')
f5 = get_fast5_file (f5file,mode='r')
reads_in_f5 = f5.get_read_ids()
number_of_reads =  len( reads_in_f5)
if number_of_reads == 1:
    print (0)
elif number_of_reads > 1:
    print (1)
