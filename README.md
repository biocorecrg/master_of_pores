[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/nanopore.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/nanopore/builds)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/179323639.svg)](https://zenodo.org/badge/latestdoi/179323639)
[![Nextflow version](https://img.shields.io/badge/Nextflow-19.10.0-brightgreen)](https://www.nextflow.io/)
[![Singularity version](https://img.shields.io/badge/Singularity-v2.6.1-green.svg)](https://www.sylabs.io/)
[![Docker version](https://img.shields.io/badge/Docker-v19.03-blue)](https://www.docker.com/)


# ![Nanopore analysis pipeline](https://github.com/biocorecrg/nanopore_analysis/blob/master/docs/logo_master.jpg) 

<img align="right" href="https://biocore.crg.eu/" src="https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png" />

<br/>


# Nanopore analysis pipeline
Nextflow pipeline for analysis of Nanopore reads (from RNA/cDNA/DNA). This project is in collaboration with [Eva Novoa's group](https://www.crg.eu/en/programmes-groups/novoa-lab). If you use this tool you can cite our pre-print:

[Parallel and scalable workflow for the analysis of Oxford Nanopore direct RNA sequencing datasets
Luca Cozzuto, Huanle Liu, Leszek P. Pryszcz, Toni Hermoso Pulido, Julia Ponomarenko, Eva Maria Novoa
doi: https://doi.org/10.1101/818336](https://www.biorxiv.org/content/10.1101/818336v1)


## Pre-requisites
For using the pipeline [Nextflow](https://www.nextflow.io/) and a linux container engine (either [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/guides/3.1/user-guide/cli/singularity_apps.html)) need to be installed. 
The pipeline can be run in Mac OSX and Linux operative systems.  

### Installation
For installing Nextflow:

```bash
curl -s https://get.nextflow.io | bash
```

The pipeline can be cloned in this way using **git**:

```bash
git clone --depth 1 https://github.com/biocorecrg/master_of_pores.git
```

Because of redistribution restriction of the basecallers **Albacore** and **Guppy** we cannot provide them inside the docker image, so you would need to download the binaries from the official website https://nanoporetech.com and place them inside the **master_of_pores/bin** folder.

#### Albacore
Download the whel file.

```bash
pip3 install --target=./albacore ont_albacore-2.1.7-cp36-cp36m-manylinux1_x86_64.whl
$ ln -s albacore/bin/multi_to_single_fast5 
$ ln -s albacore/bin/read_fast5_basecaller.py
```
#### Guppy
There are two version fo Guppy, one that runs on CPUs and the other works on both CPUs and GPUs. The difference of speed between GPUs and GPU is more than 10 times.

```bash
cd master_of_pores/bin
tar -zvxf ont-guppy_3.1.5_linux64.tar.gz
ln -s ont-guppy_3.1.5_linux64/ont-guppy/bin/guppy_* .
````
In case you want to use the GPU you need to install the CUDA drivers:
https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html 

### Running the pipeline

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

```nextflow run main.nf -with-docker
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

-----
## Pipeline workflow

<img align="middle" src="https://raw.githubusercontent.com/biocorecrg/master_of_pores/master/docs/dag_graph.png" />
