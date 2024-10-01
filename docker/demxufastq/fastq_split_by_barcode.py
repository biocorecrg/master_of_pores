#!/usr/bin/env python3
desc="""Split FastQ by barcode
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Modified by Luca Cozzuto for MoP3
"""

import gzip, os, sys
import pandas as pd
from datetime import datetime
from Bio import SeqIO

def split_by_barcode(demux, fastq, outname, minbaseQ):
	"""Split FastQ file by barcode"""
	outdir = os.path.dirname(outname)
	if outdir and not os.path.isdir(outdir):
		os.makedirs(outdir)

	df = pd.read_csv(demux, sep="\t")
	read2bc = {r: b for r, b in df.loc[df["baseQ"]>=minbaseQ, ["read_id", "barcode"]].to_numpy()}

	outs = {}
	b2c = {}
	i = k = 0
	info = "{:,} processed reads; saved {:,} reads with {} barcodes: {}\r"
	for fn in fastq:
		handle = gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn, 'rt')
		for i, r in enumerate(SeqIO.parse(handle, 'fastq'), i+1):
			read_id = r.id
			if not i%10000: sys.stderr.write(info.format(i, k, len(b2c), b2c))
			if read_id not in read2bc: continue
			bc = read2bc[read_id]
			if(str(bc) != 'nan'):
				if bc not in outs:
					outs[bc] = gzip.open("%s.bc_%d.fastq.gz"%(outname, bc), "wt")
					b2c[bc] = 0
				outs[bc].write(r.format('fastq'))
				b2c[bc] += 1
				k += 1
	sys.stderr.write(info.replace('\r','\n').format(i, k, len(b2c), b2c))
	for bc, out in outs.items():
		out.close()

def main():
	import argparse
	usage	= "%(prog)s -v" #usage=usage,
	parser	= argparse.ArgumentParser(description=desc, epilog=epilog, \
									  formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('--version', action='version', version='1.0a')
	parser.add_argument("-v", "--verbose", action="store_true", help="verbose")
	parser.add_argument("-i", "--demux", required=True, help="demux file name")
	parser.add_argument("-f", "--fastq", nargs="+", help="input FastQ file(s)")
	parser.add_argument("-o", "--outname", required=True,
						help="output basename (.bc_?.fq.gz will be added)")
	parser.add_argument("-b", "--minbaseQ", default=50, type=int,
						help="minimum demux quality [%(default)s]")

	o = parser.parse_args()
	if o.verbose:
		sys.stderr.write("Options: %s\n"%str(o))

	split_by_barcode(o.demux, o.fastq, o.outname, o.minbaseQ)

if __name__=='__main__':
	t0 = datetime.now()
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("\nCtrl-C pressed!		 \n")
	#except IOError as e:
	#	 sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
	dt = datetime.now()-t0
	sys.stderr.write("#Time elapsed: %s\n"%dt)
