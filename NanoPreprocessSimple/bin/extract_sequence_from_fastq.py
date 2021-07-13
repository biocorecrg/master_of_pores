#!/usr/bin/env python
import sys
import gzip
import pprint
import os

usage = '''
created by Huanle and Luca for Master of Pores! :)
python extract_sequence_from_fastq.py <table produced by demultiplexoon>  <fastq_file>
'''

if len (sys.argv) < 3:
    print (usage)
    exit(1)

def fopen (f):
    if f.endswith('.gz'):
        return (gzip.open(f,'rt'))
    else:
        return (open (f,'rt'))


IDs = dict ()
if True:
    fh = fopen (sys.argv[1])
    for l in fh:
        ary = l.strip().split()
        IDs[ary[1]] = ary[2]
fh.close()

outprefix = os.path.splitext(sys.argv[2])[0]
outext = os.path.splitext(sys.argv[2])[1]
if (outext == '.gz'): 
	outprefix = os.path.splitext(outprefix)[0]

fh = fopen (sys.argv[2])
for l in fh:
    rd = l.strip().split()[0]
    rd = rd.replace('@','')
    seq = next(fh).strip()
    inf = next (fh).strip()
    q = next(fh).strip()
    if rd in IDs:
        outfile = outprefix + "." + IDs[rd] + ".fastq"
        fw = open(outfile,"a+")
        rd = '@'+rd
        string = "{0}\n{1}\n{2}\n{3}\n".format(rd,seq,inf,q)
        fw.write(string)
        #print (string)
fh.close()
