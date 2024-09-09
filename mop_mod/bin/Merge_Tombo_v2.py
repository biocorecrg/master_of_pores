#!/usr/bin/env python

##Script to merge all Tombo results: statistics, ivt coverage, sample coverage and kmer score

#Import libraries:
import sys

import pyBigWig

#Import input files:
statistic = pyBigWig.open(sys.argv[1])
cov_ivt = pyBigWig.open(sys.argv[2]) 
cov_sample = pyBigWig.open(sys.argv[3])
output_name = sys.argv[4]

#Create output with header:
f = open("".join([output_name, '_Tombo_Output.tsv']), "w")
print('\t'.join(["Ref_Position","Chr","Position","Tombo_SiteScore","Coverage_Sample","Coverage_IVT","Tombo_KmerScore"]), file=f)

#Process results and close output file:
for transcript in statistic.chroms().keys():

    for i in range(0, len(cov_sample.intervals(transcript))):
        position = cov_sample.intervals(transcript)[i][1]
        ref_pos = "_".join([transcript, str(position)])
        stat = statistic.values(transcript, position-1, position)[0]
        ivt = cov_ivt.values(transcript, position-1, position)[0]
        samples = cov_sample.values(transcript, position-1, position)[0]
        
        #Calculate the kmer score:
        try:
            kmer_score = sum(statistic.values(transcript, position-3, position+2))
        except:
            kmer_score = "nan"

        #Output the processed data:
        row = [ref_pos, transcript, position, stat, samples, ivt, kmer_score]
        print('\t'.join([str(x) for x in row]), file=f)

f.close()