#!/usr/bin/env python
import sys, warnings,os,glob,collections

__author__ = 'Luca.Cozzuto@crg.eu'
__version__ = '0.1'
__email__ = 'same as author'

usage = '''
        python intersect_tombo.py fasta 0.5 output.txt 
'''

if len (sys.argv) < 3:
    print (usage)
    exit(0)

pattern   = sys.argv[1]
threshold = sys.argv[2]
outfile   = sys.argv[3]

files = glob.glob('*' + pattern)
ids = collections.defaultdict(dict)
seqs= collections.defaultdict(dict)
seqID= ""
nfiles = len(files)

for file in files:
	with open(file) as f:
		for line in f:
			if (line[0] == ">"):
				line = line[1:].rstrip()
				fields = line.split(" ")
				seqID = fields[0]
				score = fields[4]
				if (score>=threshold):
					if(seqID not in ids):
						ids[seqID] = 1
					else:
						ids[seqID] = ids[seqID] + 1
			else:
				seqs[seqID] = line.rstrip()

f = open(outfile,'w')
for id in ids:
	if(ids[id] == nfiles):
		f.write(">" + id + "\n" + seqs[id] + '\n')
f.close()


