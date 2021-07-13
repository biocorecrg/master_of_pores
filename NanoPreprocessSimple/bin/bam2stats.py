#!/usr/bin/env python3
# Report stats (mapped reads and identity to reference) from samtools stats 
# for bam file(s) ignoring secondary, suplementary and qc failed alignments
#
# USAGE: bam2stats.py bam1 bam2 ... bamN

import os, subprocess, sys

def bam2stats(fn, flag=3840):
    """Get stats from samtools stats"""
    args = ["samtools", "stats", "-F%s"%flag, fn]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE)
    k2v = {}
    for l in proc.stdout: 
        l = l.decode("utf-8")
        if l.startswith('SN'): 
            ldata = l[:-1].split()#; print(ldata) 
            kv = [[]]
            for e in ldata[1:]:
                kv[-1].append(e)
                if e.endswith(':'):
                    kv[-1][-1] = kv[-1][-1][:-1]
                    kv.append([])
            k2v[" ".join(kv[0])] = kv[1]
    # convert digits to int
    for k in k2v:
        if k2v[k][0].isdigit():
            k2v[k] = int(k2v[k][0])
    # report if no reads mapped
    if not k2v['reads mapped']:
        return "No reads mapped"
    text = []
    text.append("{:,}\t{:.1f}%\t{:,}\t{:.1f}%".format(k2v['reads mapped'], 100*k2v['reads mapped']/k2v['sequences'], k2v['bases mapped (cigar)'], 100*k2v['bases mapped (cigar)']/k2v['total length']))
    text.append("{:,}\t{:,}".format(k2v['average length'], k2v['maximum length']))
    text.append("{:.2f}%".format(100-100*k2v['mismatches']/k2v['bases mapped (cigar)'], )) #"identity: %.2f%"%(100-k2v['mismatches']/k2v['bases mapped (cigar)'], ))
    return "\t".join(text)
    
for fn in sys.argv[1:]:
    if os.path.isfile(fn):
        sys.stdout.write("#File name\tMapped reads\tMap %\tBases\tBases %\tAvg read length\tMax read length\tidentity\n")
        sys.stdout.write("%s\t%s\n"%(fn, bam2stats(fn)))
    
'''
CHK     4691e107        9942d94c        cd9ffd51
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN      raw total sequences:    4000
SN      filtered sequences:     0
SN      sequences:      4000
SN      is sorted:      1
SN      1st fragments:  4000
SN      last fragments: 0
SN      reads mapped:   1440
SN      reads mapped and paired:        0       # paired-end technology bit set + both mates mapped
SN      reads unmapped: 2560
SN      reads properly paired:  0       # proper-pair bit set
SN      reads paired:   0       # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      726     # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 6801
SN      total length:   136941  # ignores clipping
SN      bases mapped:   109284  # ignores clipping
SN      bases mapped (cigar):   108908  # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     14898   # from NM fields
SN      error rate:     1.367944e-01    # mismatches / bases mapped (cigar)
SN      average length: 34
SN      maximum length: 401
SN      average quality:        3.5
SN      insert size average:    0.0
SN      insert size standard deviation: 0.0
SN      inward oriented pairs:  0
SN      outward oriented pairs: 0
SN      pairs with other orientation:   0
SN      pairs on different chromosomes: 0

'''

