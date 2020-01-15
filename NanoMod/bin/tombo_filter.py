#!/usr/bin/env python
import argparse, re
import numpy as np
from collections import Counter as cnt
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument ('-t','--tomboouputs', required=True, dest='tbs', action='append', help='knockout sample tombo prediciton results')
parser.add_argument ('-p','--percentage', nargs = '?', const=0.5, default=0.5,type = float, help='threshod to filter for positives; (0,1]; default is 0.5')
parser.add_argument ('-o','--output', required = True, help= 'output file name')
parser.add_argument ('-m','--motif', type=str, default='[AGCTUagctu]',help='motif to be kept in the final results; default will keep all motifs; "-m [AG][AG]AC[ACT]" will keep RRACH motifs; ')
parser.add_argument ('-d','--num_duplicates', required = True, type=int, default=1, help='a site has to be at least predicted in this many samples to be used to determine final modificaiton status; default is 3')

args = parser.parse_args()

def read_tombo_fasta (fasta, motif):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	info = defaultdict (list)
	for f in fasta:
		with open  (f, 'r') as fh:
			for l in fh:
				if l.startswith ('>'):
					mtf = next (fh)
					ary = l.strip().split()
					id = ary[0][1:]
					strand = id.split(':')[-1]
					#id = ':'.join (id.split(':')[:2])
					prob = float (ary[-1])
					if strand == '-':
						# Are we sure that is not already reverse complement?
						mtf = "".join(complement.get(base,base) for base in reversed(mtf))
					id = ':'.join (id.split(':')[:2])
					id = id + ":" + mtf.rstrip()
					if not re.match (motif, mtf):
						continue
					info[id].append ((mtf,prob))
	return info

def score (sample_info, threshold, duplicates):
	scores = dict()
	for s in sample_info:
		num = len (sample_info[s])
		if num < duplicates:
			continue
		probs = [x[1] for x in sample_info[s]]
		cond = all ([x>threshold for x in probs])
		if (cond):
			scores[s] = 1
		else:
			scores[s] = np.mean (probs)
	return scores


if __name__ == '__main__':
	tbs = args.tbs
	motif = args.motif
	out = open (args.output,'w')
	tombo_info = read_tombo_fasta (tbs,args.motif) ; print (tombo_info)
	tombo_score = score (tombo_info, args.percentage, args.num_duplicates); print (tombo_score)
	mod = ''
	out.write ('#chr,position,motif,modification_status\n')
	for s in tombo_score:
		if tombo_score[s] > args.percentage:
			mod = 'YES'
		else:
			mod = 'NO'
		s = s.replace (':',',')
		out.write ("{},{}\n".format(s,mod))
	out.close()
