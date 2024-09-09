#!/usr/bin/env python3

desc = "Concatenate file(s) containing median current intensity values per position values."

"""
# example command:

./merge_nanopolish_duckdb.py -t 8 -i mod_batch_0.perpos_median_dset/contig_group=g1/*.parquet -o result.gz -l contig_groups.tsv -g g1

"""


# Import the required libraries:

import argparse
import os

import duckdb


def concat_csv_gz(list_of_files, output_fn):
    file_names = " ".join(list_of_files)
    cat_command = f"cat {file_names} > {output_fn}"
    os.system(cat_command)
    #remove the intermediates
    purge_command = f"rm  {file_names}"
    os.system(purge_command)


def process_contig_group(
    parquet_inputs_fn_list, num_threads, countig_groups_fn, group_id
):
    def parse_contig_groups(tsv_fn):
        my_group_contigs = []
        with open(tsv_fn) as f:
            for line in f:
                line = line.strip()
                sl = line.split()
                contig = sl[0]
                my_group_contigs.append(contig)

        assert len(my_group_contigs) > 0
        return my_group_contigs

    contig_list = parse_contig_groups(countig_groups_fn)
    contig_list_sqlish = ", ".join(f"'{w}'" for w in contig_list)

    # print(contig_list_sqlish)
    # set up DuckDB
    duckdb_pragma_1 = "PRAGMA enable_object_cache"
    duckdb_pragma_2 = f"PRAGMA threads={num_threads}"

    con = duckdb.connect(database=":memory:")
    con.execute(duckdb_pragma_1)
    con.execute(duckdb_pragma_2)

    # hardwired location of file for testing only
    # group_parquet_glob = f"*_dset/contig_group={group_id}/*parquet"
    # print(group_parquet_glob)

    # sql_option_1 = f"CREATE VIEW contig_grp AS SELECT * FROM parquet_scan({parquet_inputs_fn_list})"

    sql_option_2 = f"CREATE TABLE contig_grp AS SELECT * FROM parquet_scan({parquet_inputs_fn_list}) WHERE contig IN ({contig_list_sqlish}) ORDER BY contig, position, reference_kmer"

    if debug_me == True:
        print(option_1)
        print(option_1)
    con.execute(f"{sql_option_2}")

    # optional, needs to be benchmarked
    # sql_sort = ""

    sql_index = (
        "CREATE INDEX con_pos_kmer ON contig_grp (contig, position, reference_kmer)"
    )
    con.execute(f"{sql_index}")

    csv_gz_fn_list = []

    for contig_id in contig_list:
        out_csv_fn = f"remove-me_{contig_id}_merged.csv"
        print(out_csv_fn)
        contig_select_sql = f"COPY (SELECT contig, position, reference_kmer, median(median) AS 'median', sum(coverage) AS 'coverage' \
            FROM contig_grp \
            WHERE contig='{contig_id}' \
            GROUP BY contig, position, reference_kmer) TO '{out_csv_fn}' WITH (HEADER 1, DELIMITER ',')"
        print(contig_select_sql)
        con.execute(f"{contig_select_sql}")
        command_pigz = f"pigz -p{num_threads} {out_csv_fn}"

        os.system(f"{command_pigz}")
        csv_gz_fn_list.append(f"{out_csv_fn}.gz")
        print(csv_gz_fn_list)

    # print("debug 2", csv_gz_fn_list )

    return csv_gz_fn_list


def main():
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        "-i", "--input", nargs="+", help="list of parquet fn to process"
    )
    parser.add_argument("-o", "--output", help="output csv.gz")
    parser.add_argument("-t", "--threads", help="threads used by DuckDB")
    parser.add_argument(
        "-l", "--countig_groups", help="TSV file contig_name contig_group"
    )
    parser.add_argument("-g", "--group_id", help="contig_group id to process")

    a = parser.parse_args()
    # print(a)
    # Read, parse and merge data from individual eventalign files:
    csv_files_4_concat = process_contig_group(
        a.input, a.threads, a.countig_groups, a.group_id
    )
    # print("debug 1", csv_files_4_concat)
    concat_csv_gz(csv_files_4_concat, a.output)


if __name__ == "__main__":
    debug_me = False
    main()
