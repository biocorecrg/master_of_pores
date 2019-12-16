---
layout: page
title: Home
navigation: 2
---

## Pre-requisites
For using the pipeline [Nextflow](https://www.nextflow.io/) and a linux container engine (either [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/guides/3.1/user-guide/cli/singularity_apps.html)) need to be installed. 
The pipeline can be run in Mac OSX and Linux operative systems.  

## Installation


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
