---
layout: page
title: Installation
navigation: 2
---

## Pre-requisites
For using the pipeline [Nextflow](https://www.nextflow.io/) and a linux container engine (either [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/guides/3.1/user-guide/cli/singularity_apps.html)) need to be installed. 
The pipeline can be run in Mac OSX and Linux operative systems.  

## Installation


### 1. Install Nextflow (version 19.10.0)

```bash
curl -s https://get.nextflow.io | bash
```

### 2. Clone the MasterOfPores repository
The pipeline can be cloned in this way using **git**:

```bash
git clone --depth 1 https://github.com/biocorecrg/master_of_pores.git
```

### 3. Install Docker and/or Singularity 
- Docker: https://docs.docker.com/install/ (version 10.03 or later is required)
- Singularity: https://sylabs.io/guides/2.6/user-guide/quick_start.html#quick-installation-steps (version 2.6.1 or 3.2.1 is required)

### 4. Download Nanopore base-calling algorithms
Because of redistribution restriction of the basecallers **Albacore** and **Guppy** we cannot provide them inside the docker image, so you would need to download the binaries from the official website https://nanoporetech.com and place them inside the **master_of_pores/NanoPreprocess/bin** folder.

#### a) Both Albacore and Guppy
```bash
cd master_of_pores/NanoPreprocess/bin
tar -zvxf ont-guppy_3.1.5_linux64.tar.gz
ln -s ont-guppy_3.1.5_linux64/ont-guppy/bin/guppy_* .
pip3 install --target=./albacore ont_albacore-2.1.7-cp36-cp36m-manylinux1_x86_64.whl
ln -s albacore/bin/read_fast5_basecaller.py .
```

#### b) Albacore
Download the wheel file.

```bash
cd master_of_pores/NanoPreprocess/bin
pip3 install --target=./albacore ont_albacore-2.1.7-cp36-cp36m-manylinux1_x86_64.whl
$ ln -s albacore/bin/multi_to_single_fast5 
$ ln -s albacore/bin/read_fast5_basecaller.py .
```
#### c) Guppy
Please note Guppy versions older than 3.1 (e.g. 3.0.3) only runs on CPUs.
Newer versions (e.g. 3.1.5 and above) works on both CPUs and GPUs. The difference of speed between CPUs and GPU is more than 10 times.

```bash
cd master_of_pores/NanoPreprocess/bin
tar -zvxf ont-guppy_3.1.5_linux64.tar.gz
ln -s ont-guppy_3.1.5_linux64/ont-guppy/bin/guppy_* .
````

### 5. Optional step: install CUDA drivers (only needed for GPU support): 

In case you want to use the GPU you need to install the [CUDA drivers](
https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html) 

### 6. Run the pipeline:
Using Singularity:
```bash
cd master_of_pores/NanoPreprocess/
nextflow run preprocessing.nf -with-singularity
```
Using Docker:
```bash
cd master_of_pores/NanoPreprocess/
nextflow run preprocessing.nf -with-docker
``` 

