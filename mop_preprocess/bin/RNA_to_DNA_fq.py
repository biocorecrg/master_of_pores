#!/usr/bin/env python 
__author__ = 'luca.cozzuto@crg.eu'
# -*- coding utf-8 -*-

import gzip
import optparse
import os
import re
#MODULES
import sys


#BODY FUNTIONS
def options_arg():
	usage = "usage: %prog -i <input fastq file> -o <output fastq file>"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option('-i', '--input', help='Input vcf file', dest="input" )
	parser.add_option('-o', '--output',help='ouput vcf File', dest="wotus" )
	(opts,args) = parser.parse_args()
	if opts.input and opts.wotus:pass
	else: parser.print_help()
	return (opts)

def __main__ ():
	parsefile(opts.input, opts.wotus)

def parsefile(file, ofile):
	count = 0
	infile = open(file, 'r')
	fwrite = open(ofile, 'a+')
	if (file.endswith('.gz')):
		infile = gzip.open(file, 'rt')
	
	for line in infile:
		count = count + 1
		if (count%4==2):
			line = line.translate(str.maketrans({'U': 'T'}))
		fwrite.write(line)

#Calling
opts = options_arg()
if opts.input and opts.wotus:
	__main__()


#RNA = a.translate(str.maketrans({'U': 'T'}))
