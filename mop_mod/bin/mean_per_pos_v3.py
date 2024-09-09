#!/usr/bin/env python

"""
speeding up the 
/usr/bin/time ./mean_per_pos_v3.py -i mod_batch_0.fast5_event_align.tsv.gz -o mod_batch_0.fast5_event_align.polars --mean
Analysing data - position level - mean
Saving results to: mod_batch_0.fast5_event_align.polars_processed_perpos_mean.parquet
1.41user 1.14system 0:00.69elapsed 366%CPU (0avgtext+0avgdata 905180maxresident)k
0inputs+88outputs (6major+14731minor)pagefaults 0swaps

"""


desc= "Calculate mean current analysis per position from nanopolish event align output."

# Import required libraries:

import argparse

import pandas as pd
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq


def parse_input(input_fn):
    
    raw_data = pl.read_csv(input_fn, sep="\t", 
    columns = ["contig", "position", "reference_kmer", "read_name", "event_level_mean"], 
    dtype={"contig": pl.Utf8, 
    "position": pl.Int32, 
    "reference_kmer": pl.Utf8, 
    "read_name": pl.Utf8, 
    "event_level_mean": pl.Float32   })
   
    return raw_data


def mean_perpos (raw_data_polars, output_prefix):
    
    output_fn =  f"{output_prefix}_processed_perpos_mean.parquet"
    
    #Calculate mean per positions:
    print('Analysing data - position level - mean')
    q = (raw_data_polars.lazy()
     .groupby( ["contig", "position", "reference_kmer"])
     .agg([
         (pl.col("event_level_mean").mean().alias("mean")), 
         (pl.col("read_name").n_unique().alias("coverage") )
         ])
     .sort(["contig", "position","reference_kmer" ])
     )
    
    result_df = q.collect()
    result_df["read_name"] = pd.Series([1 for x in range(result_df.height)])
    result_df = result_df["contig", "position", "reference_kmer","read_name","mean", "coverage"]
    result_df.to_parquet(output_fn)

    #Output parquet file:
    print(f"Saving results to: {output_fn}")

    
def median_perpos (raw_data_polars, output_prefix):
    
    output_fn =  f"{output_prefix}_processed_perpos_median.parquet"

    #Calculate mean per positions:
    print('Analysing data - position level - median')
    #sliced_data['read_name'] = 1
    q = (raw_data_polars.lazy()
     .groupby( ["contig", "position", "reference_kmer"])
     .agg([
         (pl.col("event_level_mean").median().alias("median")), 
         (pl.col("read_name").n_unique().alias("coverage") )])
     .sort(["contig", "position","reference_kmer" ])
    )
    
    result_df = q.collect()
    result_df["read_name"] = pd.Series([1 for x in range(result_df.height)])
    result_df = result_df["contig", "position", "reference_kmer","read_name","median", "coverage"]
    result_df.to_parquet(output_fn)
    #q.collect().to_parquet(output_fn, compression='zstd')
    
    #Output .csv files:
    print(f"Saving median_perpos results to: {output_fn}")
    

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
    
    #parser.add_argument("-s", "--chunk_size", default=100000, type=int, help='Size for input subsetting [%(default)s]')
    parser.add_argument("--read_level", action='store_true', help='Analysis at per read level')
    parser.add_argument("--mean", action='store_true', help='Analysis using the mean instead of the median.')

    a = parser.parse_args()

    #Process input:
    pd.set_option('display.precision', 2)
    raw_import = parse_input(a.input)

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
