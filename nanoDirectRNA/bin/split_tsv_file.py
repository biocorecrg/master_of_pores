#!/usr/bin/env python3
import os, sys, gzip, bz2

'''
Usage:
    python3 <split_bam_tsv_to_small_files.py> big.tsv number-of-reads-in-eac-small-file
    by defautlt 20000 reads will be allocated into each small file
    let huanle (elzedliu@gmail.com) know if you run into any type of errors
'''
def openfile(f):
    if f.endswith ('.gz'):
        fh = gzip.open (f,'rt')
    elif f.endswith ('bz') or f.endswith ('bz2'):
        fh = bz2.open(f,'rt')
    else:
        fh = open(f,'rt')
    return fh

num_of_reads_in_each_file = 6000
if len (sys.argv) < 2:
    print (usage)
    exit(1)

if os.path.isfile(sys.argv[1]):
    bam_tsv = sys.argv[1]
    out_prefix = bam_tsv.split('.')[0]
else:
    print ('bam_tsv file '+ sys.argv[1] + ' does not exist')
    exit(1)

if len (sys.argv) == 3:
    try:
        num_of_reads_in_each_file = int (sys.argv[2])
    except:
        print ('please input an integer number of reads in each small file')
        raise
smallfile = None
reads_cnt = 0
#READ_NAME      FLAG    CHROM   READ_POS        BASE    QUAL    REF_POS REF     OP
reads = set()
with openfile (sys.argv[1]) as fh:
    for l in fh:
        if l.startswith('#'):
            continue
        rd = l.split()[0]
        reads.add(rd)
        if reads_cnt % num_of_reads_in_each_file == 0:
            if smallfile:
                smallfile.close()
            small_filename =  out_prefix+ '_small_{}.tsv'.format(reads_cnt + num_of_reads_in_each_file)
            smallfile = open (small_filename,'w')
        smallfile.write (l)
        reads_cnt = len(reads) + 1
