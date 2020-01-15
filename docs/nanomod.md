---
layout: page
title: NanoMod
navigation: 5
---

# NanoMod
This module allows to predict the loci with RNA modifications. Data produced by NanoPreprocess are needed and in particular the reads must be aligned to the transcriptome.

## Workflow

<img src="https://raw.githubusercontent.com/biocorecrg/master_of_pores/master/docs/dag_mod.png" width="600" align="middle">


* **index_reference** index the reference file for Epinano
* **call_variants** uses Samtools for calling the variants for Epinano
* **calc_var_frequencies** it uses TSV_to_Variants_Freq.py3 for calculating the frequencies of each variants for Epinano
* **predict_with_EPInano** It predicts the modifications with Epinano
* **filter_EPInano_pred** It filers the results from Epinano using replicates if avialable

* **resquiggling** resquiggle fast5 files for Tombo
* **getModifications** it estimates the modifications using Tombo comparing WT vs KO
* **cross_tombo_pred** it gets the intersection between differen replicates
* **join_results** it gets the intersection between the Epinano and Tombo predictions

## Input Parameters
1. **input_folders** path to the folders produced by NanoPreprocessing step.
1. **comparison** tab separated text file containing the list of comparison. Here an example:
```bash
WT1 KO1
WT2 KO2
WT3 KO3
```
1. **reference** reference transcriptome
1. **output** folder
1. **coverage** read coverage threshold for prediction
1. **tombo_opt** options for tombo
1. **epinano_opt** options for epinano
1. **tombo_score** score for filtering reliable modifications (from 0.5 to 1)
1. **epinano_score** coverage score for epinano for filtering reliable modifications (integer)
1. **mofit** motif to be found (example: "[AG][AG]AC[ACT]")
1. **wt_num** number of wt samples containing that modification 
1. **ko_num** number of ko samples containing that modification (only Epinano)
1. **email**

## Results
Three folders are produced by this module:

1. Epinano, containing the results obtained with this method. You have a single file with already filtered modifications. 

```bash
geneA,126771,GAACT,5.0,0.7865299392210111,6.0,7.650999662942007e-06,YES
geneA,139467,AGACA,26.0,1.2631007354786662e-05,34.0,9.22202894828834e-14,YES
geneA,139625,AGACA,17.0,0.012299404049052885,20.0,4.64573431891912e-06,YES
geneA,192033,AGACC,11.0,1.849874054901369e-12,11.0,3.00000089999998e-14,YES
geneA,192201,AGACA,14.0,0.01469732206992497,16.0,3.00000089999998e-14,YES
...
```
2. Tombo, containing the results obtained with this method. You have one file for each comparison WT vs KO and a final one, **tombo_all.txt**, with the intersection after filtering per score. Here an example of tombo_all.txt file:

```bash
>geneA:549289:+
CTGAC
>geneA:478105:+
GAGCT
>geneA:426607:-
TTTTT
...
```

3. Comb_mod, containing the all the modifications found in Epinano and Tombo called **RNA_modifications.txt** and a Venn Diagram. Here an example of RNA_modifications.txt file:


```bash
"positions"	"epinano"	"tombo"
"chrIX-290356"	"1"	"1"
"chrX-274861"	"1"	"1"
"chrXV-513509"	"1"	"1"
"chrI-126771"	"1"	"0"
```
