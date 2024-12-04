# MoP4 - Master of Pores 4
[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/nanopore.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/nanopore/builds)
[![mop2-CI](https://github.com/biocorecrg/MOP4/actions/workflows/build.yml/badge.svg)](https://github.com/biocorecrg/MOP4/actions/workflows/build.yml)[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow version](https://img.shields.io/badge/Nextflow-24.04.4-brightgreen)](https://www.nextflow.io/)
[![Nextflow DSL2](https://img.shields.io/badge/Nextflow-DSL2-brightgreen)](https://www.nextflow.io/)
[![Singularity version](https://img.shields.io/badge/Singularity-v3.2.1-green.svg)](https://www.sylabs.io/)
[![Docker version](https://img.shields.io/badge/Docker-v20.10.8-blue)](https://www.docker.com/)

<br/>


<img align="right" href="https://biocore.crg.eu/" src="https://raw.githubusercontent.com/CRG-CNAG/BioCoreMiscOpen/master/logo/biocore-logo_small.png" />


Master of Pores is a pipeline written in Nextflow DSL2 to analyze Nanopore data.
<br/>

It can handle reads from direct RNAseq, cDNAseq, DNAseq etc.

<br/>


![MOP4](https://github.com/biocorecrg/MoP4/blob/master/img/ssj4.png?raw=true)

The name is inspired by the Metallica's [Master Of Puppets](https://www.youtube.com/watch?v=S7blkui3nQc), the logo and the images are an homage to [Akira Toriyama](https://en.wikipedia.org/wiki/Akira_Toriyama)

## Install
Please install nextflow and singularity or docker before.

Then download the repo:

```
git clone --recurse-submodules https://github.com/biocorecrg/master_of_pores.git
```


## Documentation
The documentation is available at https://biocorecrg.github.io/master_of_pores/MOP4-dev/

## Contact
Please open an issue if you encounter any issues / troubles.
However, please go over the previous issues (including closed issues) before opening a new issue, as your same exact question might have been already answered previously. Thank you!


## Reference
If you use this tool, please cite our papers:

["Nanopore Direct RNA Sequencing Data Processing and Analysis Using MasterOfPores"
Cozzuto L, Delgado-Tejedor A, Hermoso Pulido T, Novoa EM, Ponomarenko J. *N. Methods Mol Biol. 2023*;2624:185-205. doi: 10.1007/978-1-0716-2962-8_13.](https://link.springer.com/protocol/10.1007/978-1-0716-2962-8_13)

["MasterOfPores: A Workflow for the Analysis of Oxford Nanopore Direct RNA Sequencing Datasets"
Luca Cozzuto, Huanle Liu, Leszek P. Pryszcz, Toni Hermoso Pulido, Anna Delgado-Tejedor, Julia Ponomarenko, Eva Maria Novoa.
*Front. Genet., 17 March 2020.* https://doi.org/10.3389/fgene.2020.00211](https://www.frontiersin.org/articles/10.3389/fgene.2020.00211/full)
