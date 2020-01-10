---
layout: page
title: NanoTail
navigation: 4
---

## NanoTail
This module allows to estimates polyA sizes by using two different methods. Data produced by NanoPreprocess are needed and in particular the read counts / assignment must be given.

### Steps

 1. **check_reference** It verifies whether the reference is zipped and eventually unzip it
 1. **tailfindr** it runs *tailfindr* tool in parallel.
 1. **collect_tailfindr_results** It collects the results of tailfindr.
 1. **filter_bam** Bam files are filtered with *samtools* to keep only mapped reads and remove secondary alignments
 1. **tail_nanopolish** It runs *nanopolish* in parallel.
 1. **collect_nanopolish_results** It collects the results of tail_nanopolish. 
 1. **join_results** It merges the results from the two algorithms and make a plot of the correlation.


### Input Parameters

1. **input_folders** path to the folders produced by NanoPreprocessing step.
1. **nanopolish_opt** options for the nanopolish program
1. **tailfindr_opt** options for the tailfindr program
1. **reference** reference genome / transcriptome
1. **output** folder
1. **email** 


### Results
Three folders are created by the pipeline within the output folder:
1. NanoPolish: contains the output of *nanopolish* tool.
1. Tailfindr: contains the output of *tailfindr* tool.
1. PolyA_final: contains the txt files with the combined results (i.e. predicted polyA sizes). Here an example of a test:

```bash
"Read name"	"Tailfindr"	"Nanopolish"	"Gene Name"
"013a5dde-9c52-4de1-83eb-db70fb2cd130"	52.16	49.39	"YKR072C"
"01119f62-ca68-458d-aa1f-cf8c8c04cd3b"	231.64	274.28	"YDR133C"
"0154ce9c-fe6b-4ebc-bbb1-517fdc524207"	24.05	24.24	"YFL044C"
"020cde28-970d-4710-90a5-977e4b4bbc46"	41.27	56.79	"YGL238W"
```
A plot is also produced for showing the correlation between the two methods.
