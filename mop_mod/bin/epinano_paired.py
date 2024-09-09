#!/usr/bin/env python
# -- coding: utf-8 -

import argparse
import gzip
import re
from collections import Counter as cnt
from collections import defaultdict

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument ('-k','--knockout', required=True, dest='kos', action='append', help='knockout sample  epinano prediciton results')
parser.add_argument ('-w','--wildtype', required=True, dest='wts', action='append', help='wildtype sample epinano prediction results')
parser.add_argument ('-c','--coverage', nargs = '?', const=5, default=5,type = int, help='minimum coverage to be considered as valid, default is 5')
parser.add_argument ('-o','--output', required = True, help= 'output file name')
parser.add_argument ('-m','--motif', type=str, default='[AGCTUagctu]',help='motif to be kept in the final results; default will keep all motifs; "-m [AG][AG]AC[ACT]" will keep RRACH motifs; ')
parser.add_argument ('-dk','--ko_duplicates', required = True, type=int, default=1, help='a site has to be at least predicted in this many samples to be used to determine final modificaiton status; default is 1')
parser.add_argument ('-dw','--wt_duplicates', required = True, type=int, default=1, help='a site has to be at least predicted in this many samples to be used to determine final modificaiton status; default is 1')
parser.add_argument ('-gz', '--gzipped', default=False, dest='isgzip',  action='store_true', help='indicates if the input files are gzipped or not')
args = parser.parse_args()

'''
python3 compare_wt_ko_predicitons.py -k ko.csv -k ko3.csv -k ko2.csv -w wt.csv -c 5 -o out -m [AG][AG]AC[ACT]
'''

def read_result (cov, mtf, epinano_result_files, gzip):
	'''
	input1: SVM.py predicted results
	input2: minimum depth
	input3: motif to be saved; default is RRACH
	'''
	motif = dict()
	info = defaultdict(list)
	for f in epinano_result_files:
		with fopen (f, gzip) as fh:
			for l in fh:
				if l.startswith('#'):
					continue
				ary = l.strip().split(',')
				mapping_cov = float(ary[3].split(':')[2])
				if mapping_cov < cov:
					continue
				if not re.match (mtf,ary[0]):
					continue
				win,ref,kmer,probm = ary[1],ary[2],ary[0],float(ary[-2])
				k = ref + ',' + win.split(':')[2]
				info[k].append ((mapping_cov,probm))
				motif[k] = kmer
	return info, motif

def fopen (file, isgzip):
	fh = ""
	if isgzip:
		fh = gzip.open(file, 'rt', encoding='utf-8')
	else:
		fh = open (file)
	return fh


def scoring (common_sites, sample_info, repeats):
	'''
	score prediction results
	if it is one on one comparison:
		use the probm as it is
	else:
		if (s1 ≥ 0.5 and s2 ≥ 0.5 and .. sn ≥ 0.5):
			M = 1
		else:
			M = (s1 + s2 + sn)/n
		if (Mwt/Mko) > 1.5 and Mwt > 0.5:
			status = modified
		else:
			status = unmodified
	'''
	scores = dict ()
	probs = defaultdict(list)
	for s in common_sites:
		probms = [x[1] for x in sample_info[s]]
		if len (probms) < repeats:
			continue
		elif len (probms) == 1:
			scores[s] = np.mean (probms)
		else:
			cond = [x > 0.5 for x in probms]
			if all (cond):
				scores[s] = 1
			else:
				scores[s] = np.mean (probms)
		probs[s] = probms
	return scores,probs

if __name__ == '__main__':
	kos = args.kos
	wts = args.wts
	cov = args.coverage
	ko_duplicates = args.ko_duplicates
	wt_duplicates = args.wt_duplicates
	isgzip = args.isgzip
	motif = args.motif
	out = open (args.output, 'w')

	head = "#chr,position,motif,modification_status\n"
	out.write(head)
	wt_inf, mtfs = read_result (cov,motif,wts,isgzip)
	ko_inf, mtfs = read_result (cov,motif,kos,isgzip)
	all_wt_sites_lst = [x for x in wt_inf]
	all_ko_sites_lst = [x for x in ko_inf]
	commonSites = set (all_wt_sites_lst).intersection (set(all_ko_sites_lst))
	wt_scores, wt_probs = scoring (commonSites, wt_inf, wt_duplicates)
	ko_scores, ko_probs = scoring (commonSites, ko_inf, ko_duplicates)
	for s in commonSites:
		if not (s in wt_scores and s in ko_scores):
			continue
		if wt_scores[s]/ko_scores[s] > 1.5 and wt_scores[s] > 0.5:
			out.write ('{},{},{}\n'.format(s,mtfs[s],"YES"))
		else:
			out.write ('{},{},{}\n'.format(s,mtfs[s],"NO"))
	out.close()
