.. _home-page-about:

*******************
About Master of Pores
*******************

.. autosummary::
   :toctree: generated

Master of Pores is a pipeline writte in Nextflow DSL2 for the analysis of Nanopore data. It can handle reads from direct RNAseq, cDNAseq, DNAseq etc.

The pipeline is composed by four modules:
   - mop_preprocess: preprocessing
   - mop_mod: detecting chemical modifications. It reads the output directly from mop_preprocess
   - mop_tail: estimating polyA tail size. It reads the output directly from mop_preprocess
   - mop_consensus: it generates a consensus from the predictions from mop_mod. It reads the output directly from mop_mod


The name is inspired by Metallica's `Master Of Puppets <https://www.youtube.com/watch?v=S7blkui3nQc>`_

.. image:: ../img/master_red.jpg
  :width: 400  

This is a joint project between `CRG bioinformatics core <https://biocore.crg.eu/>`_ and `Epitranscriptomics and RNA Dynamics research group <https://public-docs.crg.es/enovoa/public/website/index.html>`_.


Reference
======================

If you use this tool, please cite our paper:

`"MasterOfPores: A Workflow for the Analysis of Oxford Nanopore Direct RNA Sequencing Datasets" <https://www.frontiersin.org/articles/10.3389/fgene.2020.00211/full>`_ Luca Cozzuto, Huanle Liu, Leszek P. Pryszcz, Toni Hermoso Pulido, Anna Delgado-Tejedor, Julia Ponomarenko, Eva Maria Novoa. Front. Genet., 17 March 2020.



