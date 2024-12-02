.. _home-page-mopconsensus:

*******************
MOP_CONSENSUS
*******************

.. autosummary::
   :toctree: generated

This pipeline takes as input the output from MOP_MOD with all the four worklows. It outputs the consensus of the diferent predictions running the tool `Nanoconsensus <https://github.com/ADelgadoT/NanoConsensus>`__ in parallel on each transcript for each comparison.


Input Parameters
======================

The input parameters are stored in yaml files like the one represented here:

.. literalinclude:: ../mop_consensus/params.yaml
   :language: yaml


How to run the pipeline
=============================

Before launching the pipeline,user should:

1. Decide which containers to use - either docker or singularity **[-with-docker / -with-singularity]**.
2. Fill in both **params.config** and **tools_opt.tsv** files.

To launch the pipeline, please use the following command:

.. code-block:: console

   nextflow run mop_consensus.nf -params-file params.yaml  -with-singularity > log.txt


You can run the pipeline in the background adding the nextflow parameter **-bg**:

.. code-block:: console

   nextflow run mop_consensus.nf -params-file params.yaml -with-singularity -bg > log.txt

You can change the parameters either by changing **params.config** file or by feeding the parameters via command line:

.. code-block:: console

   nextflow run mop_consensus.nf -params-file params.yaml -with-singularity -bg --output test2 > log.txt


You can specify a different working directory with temporary files:

.. code-block:: console

   nextflow run mop_consensus.nf -params-file params.yaml -with-singularity -bg -w /path/working_directory > log.txt


Results
====================

Here an example of a result:

.. image:: ../img/nanocons.png
  :width: 800
