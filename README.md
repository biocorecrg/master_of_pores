[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/nanopore.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/nanopore/builds)
[![Singularity version](https://img.shields.io/badge/Singularity-v2.6.1-green.svg)](https://www.sylabs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.org/biocorecrg/master_of_pores.svg?branch=master)](https://travis-ci.org/biocorecrg/master_of_pores)
[![DOI](https://zenodo.org/badge/179323639.svg)](https://zenodo.org/badge/latestdoi/179323639)


# ![Nanopore analysis pipeline](https://github.com/biocorecrg/nanopore_analysis/blob/master/docs/logo_master.jpg) 

<img align="right" href="https://biocore.crg.eu/" src="https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png" />

<br/>


# Nanopore analysis pipeline
Nextflow pipeline for analysis of Nanopore reads (from RNA/cDNA/DNA). This project is in collaboration with [Eva Novoa's group](https://www.crg.eu/en/programmes-groups/novoa-lab).

## Docker image
The docker images were generated starting from the Docker files in this repository and uploaded to DockerHub: [**Nanopore image**](https://hub.docker.com/r/biocorecrg/nanopore) [**Porechop image**](https://hub.docker.com/r/biocorecrg/porechop)

Docker images can be converted to singularity images using the command:

```bash
singularity pull docker://biocorecrg/nanopore:0.8
singularity pull docker://biocorecrg/porechop:0.1
```

**Singularity** version >= **2.6.1** is needed.

You can create an environmental variable:
```bash
export RUN="singularity exec --cleanenv nanopore-0.8.simg"
export PORECHOP="singularity exec --cleanenv porechop-0.1.simg"
```

and then use it like this:

```bash
$RUN python --version

Python 3.6.3
```

## The pipeline
The pipeline can be cloned in this way:

```bash
git clone https://github.com/biocorecrg/nanopore_analysis.git
```

The pipeline is composed of two modules:
- Preprocessing: 
  - preprocessing.nf 
- RNA modifications: (Future plan)
  - to be decided

### preprocessing.nf
Input files are either multifast5 or single fast5 files containing reads from genomic DNA, cDNA or RNA sequencing. 
It needs a reference sequence (genome or transcriptome).




  - baseCalling with **Albacore** or **Guppy**
  - demultiplexing (optional) with **Guppy** or **porechop** 
  - tarFast5
  - QC with **MinIONQC.R**
  - **fastQC**
  - mapping with **minimap2** or **graphmap**. **Samtools** is the used for conversion.
  - alnQC2 with custom script **bam2stats**.
  - alnQC2 with **NanoPlot**.
  - **multiQC** for the final report
  

You can launch the pipeline in this way choosing either the parameter **-with-singularity** or **with-docker** depending on which containers you want to use:

```bash

nextflow run -bg preprocessing.nf -with-singularity

N E X T F L O W  ~  version 19.01.0
Launching `preprocessing.nf` [wise_colden] - revision: 6a828b7af6


╔╦╗┌─┐┌─┐┌┬┐┌─┐┬─┐  ┌─┐┌─┐  ╔═╗╔═╗╦═╗╔═╗╔═╗
║║║├─┤└─┐ │ ├┤ ├┬┘  │ │├┤   ╠═╝║ ║╠╦╝║╣ ╚═╗
╩ ╩┴ ┴└─┘ ┴ └─┘┴└─  └─┘└    ╩  ╚═╝╩╚═╚═╝╚═╝
                                                                                       
====================================================
BIOCORE@CRG Preprocessing of Nanopore data (gDNA, cDNA or RNA) - N F  ~  version 0.1
====================================================

kit                       : SQK-RNA001
flowcell                  : FLO-MIN106
fast5                     : ./data_2/RNAAB056712_wt2/*.fast5
multi5                    : YES
reference                 : /nfs/software/bi/biocore_tools/git/nextflow/master_of_pores/anno/genome.fa.gz

seqtype                   : RNA
output                    : out_RNAAB056712_wt2
granularity               : 2

basecaller                : guppy
basecaller_opt            : 
GPU                       : OFF
demultiplexing            :  
demultiplexing_opt        :  
barcodekit                : EXP-NBD104
mapper                    : minimap2
mapper_opt                : 
map_type                  : spliced

email                     : luca.cozzuto@crg.eu
```

You can change them by editing the **preproc.config** file or using the command line (each param name needs to have the characters **--** before): 

```bash
nextflow run preprocessing.nf -with-singularity -bg --granularity 20 > log.txt
```

To resume a previous execution that failed at a certain step or if you change a parameter that affects only some steps you can use the **Netxtlow** parameter **-resume** (only one dash!):


```bash
nextflow run preprocessing.nf -with-singularity -bg -resume > log.txt

...

[warm up] executor > crg
[e8/2e64bd] Cached process > baseCalling (RNA081120181_1)
[cd/396b0b] Cached process > tarFast5 (RNA081120181_1)
[b2/21f680] Cached process > QC (RNA081120181_1)
[c8/3f5d17] Cached process > mapping (RNA081120181_1)
...

```

-----------------------------------------------------

# TODO
## First module:
* De-multiplexing
* cDNA and RNA gene counts


Currently the pipeline has the following steps:

1. **baseCalling**: it converts fast5 files into a single fastq file using either **Albacore** or **Guppy** depending on the parameter **basecaller**.
1. **tarFast5**: it saves the fast5 files in a single archive using **tar**
1. **QC**: performed using **MinIONQC.R**
1. **fastQC**: QC on fastq files.
1. **mapping**: it maps either to the transcriptome or to the genome (parameter **reftype**: **T** or **G**). Alignment is then converted into a sorted bam file using **samtools**. The mapping is performed using either **minimap2** or **graphmap** (parameter **mapping**: **minimap2** or **graphmap**). In case the input sequence is **RNA** (specified by **seqtype** parameter) the **U** are converted into **T** before the alignment.
1. **alnQC**: quality control over aligned reads
1. **alnQC2**: quality control over aligned reads




1. **multiQC**: it collects data from different QC and groups them into a single report

The pipeline accept both single fast5 reads or multi-fast5. You need to specify the format using the parameter **multi5**.

-----

