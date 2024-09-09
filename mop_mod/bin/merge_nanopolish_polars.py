#!/usr/bin/env python
desc = "Concatenate file(s) containing median current intensity values per position values."

# Import required libraries:
import argparse
import glob
import os

import pandas as pd
import polars as pl
import pyarrow.parquet as pq


def process_contigs(input_parquet_files, out_dir, num_threads):
    "no contigs groups as input yet"

    # tmp dir to output partitioned by contig prquets
    datasets_dir = "out_datasets_merger"
    os.system(f"mkdir {out_dir}")
    # python set to collect contigs/transcripts names
    contigs_set = set()

    for parquet_fn in input_parquet_files:
        # print(parquet_fn)
        q = pl.scan_parquet(
            parquet_fn,
            columns=["contig", "position", "reference_kmer", "median", "coverage"],
        )
        pl_df = q.collect()
        transcripts = pl_df["contig"].unique().to_list()
        contigs_set.update(transcripts)
        pq.write_to_dataset(
            pl_df.to_arrow(), root_path=datasets_dir, partition_cols=["contig"]
        )

    # print("got here 001")

    for contig_name in contigs_set:
        parquet_files_glob = f"{datasets_dir}/contig={contig_name}/*.parquet"
        out_csv_fn = f"{out_dir}/{contig_name}_merged_nanopol.csv"

        q = (
            pl.scan_parquet(parquet_files_glob)
            .groupby(["position", "reference_kmer"])
            .agg(
                [
                    (pl.col("median").median().alias("median")),
                    (pl.col("read_name").n_unique().alias("coverage")),
                ]
            )
            .sort(["position", "reference_kmer"])
        )
        contig_df = q.collect()

        contig_df["contig"] = pd.Series([contig_name for x in range(contig_df.height)])

        contig_df.to_csv(out_csv_fn, sep="\t")
        command_pigz = f"pigz -p{num_threads} {out_csv_fn}"
        os.system(f"{command_pigz}")


def main():
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-d", "--dir", help="Input top dir to process.")
    parser.add_argument("-o", "--output", help="Output dir")
    parser.add_argument("-t", "--threads", help="threads used by pigz")

    a = parser.parse_args()
    print(a)
    input_parquet_files = glob.glob(f"{a.dir}/*.parquet")
    print(input_parquet_files)
    # Read, parse and merge data from individual eventalign files:
    process_contigs(input_parquet_files, a.output, a.threads)


if __name__ == "__main__":
    debug_me = False
    main()
