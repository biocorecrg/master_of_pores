---
layout: page
title: Running the pipeline
navigation: 3
---

# Running the pipeline
## NanoPreprocess

Input files are either multifast5 or single fast5 files containing reads from genomic DNA, cDNA or RNA sequencing. 
They will be basecalled and eventually demultiplexed and aligned to a reference sequence (genome or transcriptome).

Steps:
  - Basecalling with **Albacore** or **Guppy**
  - Demultiplexing (optional) with **Guppy** or **Porechop** 
  - Tar of fast5 files (in case they are single sequence fast5)
  - QC with **MinIONQC.R**
  - QC with **fastQC**
  - Filtering with **NanoFilt**.
  - Mapping with **minimap2** or **graphmap**. **Samtools** is the used for conversion.
  - QC of aligned reads with custom script **bam2stats**.
  - QC of aligned reads with **NanoPlot**.
  - Final report with **multiQC**
  

You can launch the pipeline choosing either the parameter **-with-singularity** or **with-docker** depending on which containers you want to use:

```
nextflow run main.nf -with-docker
N E X T F L O W  ~  version 0.31.1
Launching `biocorecrg/master_of_pores` [pensive_boyd] - revision: fc7613225b [master]


╔╦╗┌─┐┌─┐┌┬┐┌─┐┬─┐  ┌─┐┌─┐  ╔═╗╔═╗╦═╗╔═╗╔═╗
║║║├─┤└─┐ │ ├┤ ├┬┘  │ │├┤   ╠═╝║ ║╠╦╝║╣ ╚═╗
╩ ╩┴ ┴└─┘ ┴ └─┘┴└─  └─┘└    ╩  ╚═╝╩╚═╚═╝╚═╝
                                                                                       
====================================================
BIOCORE@CRG Preprocessing of Nanopore data (gDNA, cDNA or RNA) - N F  ~  version 0.1
====================================================

kit                       : SQK-RNA001
flowcell                  : FLO-MIN106
fast5                     : ./data/multifast/*.fast5
multi5                    : YES
reference                 : ./anno/curlcake_constructs.fasta.gz

seqtype                   : RNA
output                    : output
granularity               : 1
qualityqc                 : 5

basecaller                : albacore
basecaller_opt            : 
GPU                       : OFF
demultiplexing            :  
demultiplexing_opt        :  
barcodekit                : 
filter                    : nanofilt
filter_opt                : -q 0 --headcrop 5 --tailcrop 3 --readtype 1D
mapper                    : minimap2
mapper_opt                : 
map_type                  : spliced

email                     : ""
```

You can change them by editing the **preproc.config** file or using the command line (each param name needs to have the characters **--** before): 

```bash
nextflow run main.nf -with-singularity -bg --granularity 20 > log.txt
```

To resume a previous execution that failed at a certain step or if you change a parameter that affects only some steps you can use the **Netxtlow** parameter **-resume** (only one dash!):


```bash
nextflow run main.nf -with-singularity -bg -resume > log.txt

...

[warm up] executor > crg
[e8/2e64bd] Cached process > baseCalling (RNA081120181_1)
[b2/21f680] Cached process > QC (RNA081120181_1)
[c8/3f5d17] Cached process > mapping (RNA081120181_1)
...

```

-----------------------------------------------------


Currently the pipeline has the following steps:

1. **baseCalling**: it converts fast5 files into a single fastq file using either **Albacore** or **Guppy** depending on the parameter **basecaller**.
1. **QC**: performed using **MinIONQC.R**
1. **fastQC**: QC on fastq files.
1. **filtering**: Filtering of fastq files using **NanoFilt**. It is optional.
1. **mapping**: it maps either to the transcriptome or to the genome (parameter **reftype**: **T** or **G**). Alignment is then converted into a sorted bam file using **samtools**. The mapping is performed using either **minimap2** or **graphmap** (parameter **mapping**: **minimap2** or **graphmap**). In case the input sequence is **RNA** (specified by **seqtype** parameter) the **U** are converted into **T** before the alignment.
1. **alnQC**: quality control over aligned reads
1. **alnQC2**: quality control over aligned reads
1. **multiQC**: it collects data from different QC and groups them into a single report

The pipeline accept both single fast5 reads or multi-fast5. You need to specify the format using the parameter **multi5**
The parameter **granularity** is related to the amount of input file to be analyzed in a single execution. In case you have single sequence fast5 you can use a value of 2000 or up to 4000. In case you have multi-fast5 file you can go for a value of 1 or in case you use **Guppy** with GPU support a better choice can be up to 300 per time depending on the amount of GPU-RAM available. 

## NanoTail
Data produced by NanoPreprocess are needed for this module. 

Finish writing usage

## NanoMod
Data produced by NanoPreprocess are needed for this module. 

Finish writing usage

