#!/usr/bin/env python
desc= "Concatenate file(s) containing median current intensity values per position values."

# Import required libraries:
import argparse
import csv
import statistics

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pyarrow import csv


def process_dicts(files):
    intensity = dict()
    coverage = dict()
    
    for file in files:
        #Open and read input file:
            test = pq.read_table(file).to_pandas()
            #test = pq.write_table(csv_reader, '2hpf-Rep1-Batch1.parquete')
            
            #Update both dictionaries:
            for i in zip(test['contig'],test['position'],test['reference_kmer'],test['read_name'],test['median'], test['coverage']):
                
                key = ','.join(map(str, i[0:4]))
                
                #Update intensity dict:
                if key in intensity:
                    intensity[key].append(float(i[4]))
                else:
                    intensity[key] = [float(i[4])]

                #Update coverage dict:
                if key in coverage:
                    coverage[key].append(int(i[5]))
                else:
                    coverage[key] = [int(i[5])]

            #Close input file:
            #csv_file.close()
    
    #Data processing of data stored in both dictionaries:
    for key in intensity:
        intensity[key] = statistics.median(intensity[key])
        coverage[key] = sum(coverage[key])
    
    return intensity,coverage

def generate_output(output_file, intensity, coverage):
    #Create output file:
    out_file = output_file+'.tsv'
    f = open(out_file, 'w')
    
    #Print header:
    header=['contig','position','reference_kmer','read_name','median','coverage']
    print('\t'.join(header),file=f)
    
    #Print data stored in both dictionaries:
    for key in intensity:
        splitted_key = key.split(",")
        print('\t'.join(splitted_key), "{:.3f}".format(intensity.get(key)), coverage.get(key), file=f, sep='\t')
     
    f.close()
    
def main():
    parser  = argparse.ArgumentParser(description=desc)

    parser.add_argument('-i', '--input', nargs='+' ,help='Input file(s) to process.')
    parser.add_argument('-o', '--output', help='Output filename')

    a = parser.parse_args()
    
    #Read, parse and merge data from individual eventalign files:
    intensity, coverage = process_dicts(a.input)
    
    #Generate output file with processed data:
    generate_output(a.output, intensity, coverage)

if __name__=='__main__': 
    main()