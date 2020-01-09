#!/usr/bin/env python3
import argparse, re
import numpy as np
from collections import Counter as cnt
#from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument ('-k', '--knockout', required=True, help='knockout sample  epinano prediciton results')
parser.add_argument ('-w', '--wildtype', required=True, help='wildtype sample epinano prediction results')
parser.add_argument ('-c', '--coverage', nargs = '?', const=5, default=5, type = int, help='minimum coverage to be considered as valid, default is 5')
parser.add_argument ('-o', '--output', required = True, help= 'output file name')
parser.add_argument ('-m','--motif', type=str, default='[AG][AG]AC[ACT]',help='motif to be kept in the final results; default is RRACH motif aka "[AG][AG]AC[ACT]"; if this option is not supplied, all motifs will be saved ')
args = parser.parse_args()

'''
python3 final.py -k ko.csv -w wt.csv -c 5 -o oit -m [AG][AG]AC[ACT]
'''

def read_result (epinano_result_file, cov, *args):
    '''
    input1: SVM.py predicted results
    input2: minimum depth
    input3: motif to be saved; default is RRACH
    '''
    with open (epinano_result_file) as fh:
        #info = OrderedDict()
        info = dict()
        for l in fh:
            if l.startswith('#'):
                continue
            ary = l.strip().split(',')
            real_cov = float(ary[3].split(':')[2])
            if real_cov < cov:
                continue
            if args and not re.match (args[0],ary[0]):
                continue
            win,ref,kmer,probm = ary[1],ary[2],ary[0],float (ary[-2])
            k = ref + ',' + win.split(':')[2]
            info[k] = kmer,real_cov,probm
    return info


if __name__ == '__main__':
    ko = args.knockout
    wt = args.wildtype
    cov = args.coverage
    out = open (args.output, 'w')
    motif = args.motif

    head = "#chr,position,motif,cov_wt,probm_wt,cov_ko,probm_ko,modification_status\n"
    out.write(head)
    ko_inf = read_result (ko, cov, motif)
    ko_sites = set (ko_inf)
    wt_inf = read_result (wt, cov, motif)
    wt_sites = set (wt_inf)
    commonSites = ko_sites.intersection (wt_sites)
    mod = ''
    for s in sorted (commonSites):
        info_wt, info_ko = wt_inf[s], ko_inf[s]
        if info_wt[-1] / info_ko[-1] > 1.5:
            mod = 'YES'
        else:
            mod = 'NO'
        out.write ("{},{},{},{},{},{},{}\n".format(s,info_wt[0],info_wt[1],info_wt[2], info_ko[1],info_ko[2],mod))
    out.close()
