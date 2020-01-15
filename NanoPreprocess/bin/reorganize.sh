#!/bin/env bash

if [ x"$1" == x ]; then

         "please specify an output folder"

        exit 1

fi

for i in $1/fast5_files/* ; do  mkdir -p $1/`basename $i`/fast5_files;  mv $i/* $1`basename $i`/fast5_files; done
for i in $1/fast5_files/* ; do  mkdir -p $1/`basename $i`/fastq_files;  mv $1/fastq_files/*`basename $i`.fq.gz $1`basename $i`/fastq_files; done
for i in $1/fast5_files/* ; do  mkdir -p $1/`basename $i`/counts;  mv $1/counts/*`basename $i`.count $1`basename $i`/counts; done
for i in $1/fast5_files/* ; do  mkdir -p $1/`basename $i`/assigned;  mv $1/assigned/*`basename $i`.assigned $1`basename $i`/assigned; done
for i in $1/fast5_files/* ; do  mkdir -p $1/`basename $i`/alignment;  mv $1/alignment/*`basename $i`*.bam $1`basename $i`/alignment; done
rmdir $1/alignment $1/assigned $1/counts $1/fastq_files; for i in $1/fast5_files/*; do rmdir $i;  done; rmdir $1/fast5_files


