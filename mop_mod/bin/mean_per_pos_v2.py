#!/usr/bin/env python
desc= "Calculate mean current analysis per position from nanopolish event align output."

# Import required libraries:
import argparse

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def parse_input(file, size_chunks):
    
    df_chunk = pd.read_csv(file, sep='\t', chunksize=size_chunks, compression='gzip', error_bad_lines=False)
    chunk_list = list()

    # Process each portion of input file:
    for chunk in df_chunk:  
        
        chunk_filter = chunk.iloc[:,[0,1,2,3,6]]
        chunk_filter.columns = ['contig', 'position','reference_kmer', 'read_name','event_level_mean']
        chunk_filter = chunk_filter.groupby(['contig', 'position','reference_kmer', 'read_name']).agg({'event_level_mean':'mean'})
        chunk_filter.columns = ['event_level_mean']
        chunk_filter = chunk_filter.reset_index()
        
        # Once the data filtering is done, append to list
        chunk_list.append(chunk_filter)
        print('Partition {}: Processed'.format(len(chunk_list)))

    # Generate raw_input
    raw_data = pd.concat(chunk_list)
    raw_data.columns = ['contig', 'position','reference_kmer', 'read_name','event_level_mean']
    print('Concatenating partitions')

    return raw_data


def mean_perpos (sliced_data, output):
    
    #Calculate mean per positions:
    print('Analysing data - position level - mean')
    sliced_data['read_name'] = 1
    mean_perpos = sliced_data.groupby(['contig', 'position','reference_kmer', 'read_name']).agg({'event_level_mean':['mean', 'count']})
    mean_perpos.columns = ['mean', 'coverage']
    mean_perpos = mean_perpos.reset_index()

    #Output .csv files:
    print('Saving results to: {}_processed_perpos_mean.parquete'.format(output))
    pq.write_table(pa.Table.from_pandas(mean_perpos), '{}_processed_perpos_mean.parquete'.format(output))

def median_perpos (sliced_data, output):
    
    #Calculate mean per positions:
    print('Analysing data - position level - median')
    #sliced_data['read_name'] = 1
    median_perpos = sliced_data.groupby(['contig', 'position','reference_kmer'], observed=True).agg({'event_level_mean':'median', 'read_name':'nunique'})
    median_perpos.columns = ['median', 'coverage']
    median_perpos = median_perpos.reset_index()
    median_perpos['read_name'] = 1
    median_perpos = median_perpos[['contig', 'position','reference_kmer','read_name','median', 'coverage']]

    #Output .csv files:
    print('Saving results to: {}_processed_perpos_median.parquete'.format(output))
    pq.write_table(pa.Table.from_pandas(median_perpos), '{}_processed_perpos_median.parquete'.format(output))

def mean_perpos_perread (raw_data, output):

    #Calculate mean per positions:
    print('Analysing data - read level - mean')
    mean_perpos_perread = raw_data.groupby(['contig', 'position','reference_kmer', 'read_name']).agg({'event_level_mean':'mean'})
    mean_perpos_perread.columns = ['event_level_mean']
    mean_perpos_perread = mean_perpos_perread.reset_index()

    #Output .csv files:
    print('Saving results to: {}_processed_perpos_perread_mean.parquete'.format(output))
    pq.write_table(pa.Table.from_pandas(mean_perpos_perread), '{}_processed_perpos_perread_mean.parquete'.format(output))


def median_perpos_perread (raw_data, output):

    #Calculate mean per positions:
    print('Analysing data - read level - median')
    median_perpos_perread = raw_data.groupby(['contig', 'position','reference_kmer', 'read_name']).agg({'event_level_mean':'median'})
    median_perpos_perread.columns = ['median']
    median_perpos_perread = median_perpos_perread.reset_index()

    #Output .csv files:
    print('Saving results to: {}_processed_perpos_perread_median.parquete'.format(output))
    pq.write_table(pa.Table.from_pandas(median_perpos_perread), '{}_processed_perpos_perread_median.parquete'.format(output))

def main():
    parser  = argparse.ArgumentParser(description=desc)

    parser.add_argument('-i', '--input', help='Input file to process.')
    parser.add_argument('-o', '--output', help='Output filename')
    parser.add_argument("-s", "--chunk_size", default=100000, type=int, help='Size for input subsetting [%(default)s]')
    parser.add_argument("--read_level", action='store_true', help='Analysis at per read level')
    parser.add_argument("--mean", action='store_true', help='Analysis using the mean instead of the median.')

    a = parser.parse_args()

    #Process input:
    pd.set_option('display.precision', 2)
    raw_import = parse_input(a.input, a.chunk_size)

    #Analysis:

    if a.mean and a.read_level:
        mean_perpos_perread(raw_import, a.output)

    elif not a.mean and a.read_level:
        median_perpos_perread(raw_import, a.output)

    else:
        if a.mean:
            mean_perpos(raw_import, a.output)
        else:
            median_perpos(raw_import, a.output)    


if __name__=='__main__': 
    main()
