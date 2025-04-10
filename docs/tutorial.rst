.. _tutorial:

*******************
TUTORIAL
*******************

.. autosummary::
   :toctree: generated

This basic tutorial will show you hot to install and run Master of Pores in different scenarios.
  
Installing the tool and the dependencies
======================

- Install nextflow (see `here <https://www.nextflow.io/docs/latest/install.html>`_ for the full doc), java version >= 17 is required. 
                    
.. code-block:: console

   #Installing nextflow
   curl -s https://get.nextflow.io | bash

If you are on a Linux machine, you can install either `Docker <https://www.docker.com/get-started/>`_ or `singularity / apptainer <https://apptainer.org/docs/admin/main/installation.html>`_.
If you are on a Mac you can only install `Docker <https://www.docker.com/get-started/>`_.
                    
.. tip:: 
                    
    On Linux and in particular on HPC we suggest using singularity / apptainer
                    
Let's now install Master of Pores using `git clone`
                    
.. code-block:: console

  git clone --depth 1 --recurse-submodules https://github.com/biocorecrg/master_of_pores.git

  Cloning into 'master_of_pores'...
  remote: Enumerating objects: 96, done.
  remote: Counting objects: 100% (96/96), done.
  remote: Compressing objects: 100% (87/87), done.
  remote: Total 96 (delta 12), reused 56 (delta 2), pack-reused 0 (from 0)
  Receiving objects: 100% (96/96), 10.64 MiB | 14.68 MiB/s, done.
  Resolving deltas: 100% (12/12), done.
  Submodule 'BioNextflow' (https://github.com/biocorecrg/BioNextflow) registered for path 'BioNextflow'
  Cloning into '/Users/lcozzuto/ooo/master_of_pores/BioNextflow'...
  remote: Enumerating objects: 2763, done.        
  remote: Counting objects: 100% (250/250), done.        
  remote: Compressing objects: 100% (169/169), done.        
  remote: Total 2763 (delta 150), reused 163 (delta 81), pack-reused 2513 (from 2)        
  Receiving objects: 100% (2763/2763), 107.75 MiB | 10.21 MiB/s, done.
  Resolving deltas: 100% (1774/1774), done.
  Submodule path 'BioNextflow': checked out 'c70c28508dbc44c362cc77208130b24d0dbb2e78'                  

This will download the pipeline and the required submodules.
                           
Starting from fastq
======================

The test dataset is bundled with the repository. We have two small compressed fastq samples:

.. code-block:: console

  cd master_of_pores
  ls data/fastq/
  mod.fq.gz	wt.fq.gz

To analyze them, we need to go to the mop_preprocess folder and run the pipeline. All the required parameters for running the pipeline are in a yaml file. Let's check the params.yaml
                           
.. literalinclude:: ../mop_preprocess/params.yaml
   :language: yaml

The first part is for pod5 inputs, so we can ignore it. We can check the `# Needed for fastq input` part. The path of input fastq files is already specified:
                      
.. code-block:: yaml

  # Needed for fastq input
  fastq: "${projectDir}/../data/fastq/*.fq.gz"

We then need to specify the reference sequence in FASTA format and whether this is a transcriptome or a genome. In case is a genome you need to pass also the annotation in GTF format.
Then there is a section of `Actions`. You can either specify the tool for that action or turn it off using "NO" as a value.

- filtering: modifying fastq
- mapping: aligning fastq
- counting: counting read tags
- discovery: transcriptome assembly
- cram_conv: convertion of bam to cram
- subsampling_cram: subsample the bam input for generating cram

Then a new section is for specifying the `output` folder, and if you want to receive a mail or a slack message at the end of the execution. 
You need a configured mail server for sending an email and a `Slack hook <https://api.slack.com/messaging/webhooks>`_ 

Starting from pod5 (linux local and HPC)
======================


Ciao
  
.. code-block:: console

   nextflow run mop_preprocess.nf -with-singularity -params-file params.yaml > log.txt


You can run the pipeline in the background by adding the nextflow parameter **-bg**:

.. code-block:: console

   nextflow run mop_preprocess.nf -params-file params.yaml -with-singularity -bg > log.txt

You can change the parameters either by changing the yaml config file or by feeding the parameters via command line:

.. code-block:: console

   nextflow run mop_preprocess.nf -with-singularity -params-file params.yaml -bg --output test2 > log.txt


Starting from pod5 (mac, local)
======================

The command line options for each tool used in the pipeline are stored within in the same yaml file with other parameters. The section is called **progPars**. Here is an example:

.. literalinclude:: ../mop_preprocess/params.yaml
   :language: yaml
   :emphasize-lines: 44-65

The second level indicates the processing step as **basecalling** or **demultiplexing** etc, while the third indicates the tool. Finally you have the command specific command line between quotation marks.

.. note::
   You can indicate the models to be used for basecalling with dorado or dorardo-duplex as "sup,m6A_DRACH". The pipeline will try to download before and then to perform the basecalling. In case you want a specific model version you need to indicate the base simplex model as "rna002_70bps_hac@v3,pseU". You can see `here <https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models>`_ the list of models and modifications.

.. tip::
   You don't need to specify the whole path for the models of seqtagger, just the name of the model will be enough





