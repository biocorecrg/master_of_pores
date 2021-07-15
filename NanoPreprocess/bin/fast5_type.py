#!/usr/bin/env python
import sys
from ont_fast5_api.multi_fast5 import MultiFast5File
from ont_fast5_api.fast5_info import _clean
__author__ = 'Huanle.Liu@crg.eu'
__version__ = '0.2'
__email__ = 'same as author'
usage = '''
        python fast5_type.py fast5file
        return:
            0: single read fast5
            1: multi-reads fast5
'''
if len (sys.argv) !=2:
	print (usage, file=sys.stderr)
	sys.exit()
def check_file_type(f5_file):
	try:
		return _clean(f5_file.handle.attrs['file_type'])
	except KeyError:
		if len(f5_file.handle) == 0 :
			return 1 
		if len([read for read in f5_file.handle if read.startswith('read_')]) !=0 :
			return 1
		if 'UniqueGlobalKey' in f5_file.handle:
			return 0
	raise TypeError('file can not be indetified as single- or multi- read.\n' 'File path: {}'.format(f5_file.filename))
filepath = sys.argv[1]
f5_file = MultiFast5File (filepath, mode='r') 
filetype = check_file_type (f5_file)
print (filetype)
