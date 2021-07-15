## Version 1.1
* Added a new module called NanoPreprocessSimple that starts from fastq files instead of fast5 files. It allows the analysis of multiple files at a time.
* Added support to vbz compressed fast5 https://github.com/nanoporetech/vbz_compression in NanoPreprocess, NanoMod and NanoTail
* Added plots for Epinano output in NanoMod
* Added a conversion of Tombo results in bed format in NanoMod
* Added a INSTALL.sh file for automatically retrieve guppy 3.4.5 from https://mirror.oxfordnanoportal.com/, place it in NanoPreprocess/bin and making the required links
* Added profiles for being used locally and on the CRG SGE cluster


## Version 1.0
This is the original version published in the paper [MasterOfPores: A Workflow for the Analysis of Oxford Nanopore Direct RNA Sequencing Datasets](https://www.frontiersin.org/articles/10.3389/fgene.2020.00211/full)
