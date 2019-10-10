#!/bin/bash
#$ -cwd 
#$ -q biocore-el7
#$ -pe smp 1
#$ -j y
#$ -l virtual_free=80G
#$ -l h_rt=12:00:00
docker build -t biocorecrg/nanopore:0.8 .
docker push biocorecrg/nanopore:0.8
