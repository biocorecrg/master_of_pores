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

Then a new section is for specifying the `output` folder, and if you want to receive a mail or a Slack message at the end of the execution. 
You need a configured mail server for sending an `email <https://nextflow.io/docs/latest/notifications.html#mail-configuration>`_ and a `Slack hook <https://api.slack.com/messaging/webhooks>`_ 

Finally, there is a section about command-line parameters for each tool used. 

.. code-block:: yaml

  # Program params
  ProgPars:
  basecalling:
    dorado: "sup"
    dorado-duplex: "sup"
    dorado-mod: "sup,m6A_DRACH"
  demultiplexing:
  ...

To run the pipeline, just type:

.. code-block:: console

   nextflow run mop_preprocess.nf -with-docker -params-file params.yaml

You will get this as output. 

.. code-block:: console
   
    N E X T F L O W   ~  version 25.02.3-edge
   
   Launching `mop_preprocess.nf` [loving_hugle] DSL2 - revision: 3bc7696a53
   
   
   
   ====================================================
   ╔╦╗╔═╗╔═╗  ╔═╗┬─┐┌─┐┌─┐┬─┐┌─┐┌─┐┌─┐┌─┐┌─┐
   ║║║║ ║╠═╝  ╠═╝├┬┘├┤ ├─┘├┬┘│ ││  ├┤ └─┐└─┐
   ╩ ╩╚═╝╩    ╩  ┴└─└─┘┴  ┴└─└─┘└─┘└─┘└─┘└─┘
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣷⣶⣤⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⣿⣿⣿⣿⣷⡒⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⣿⣿⣿⣿⣿⣆⠙⡄⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣤⣤⣤⣤⣤⣤⣤⣤⣤⠤⢄⡀⠀⠀⣿⣿⣿⣿⣿⣿⡆⠘⡄⠀⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⢿⣿⣿⣿⣿⣿⣿⣿⣦⡈⠒⢄⢸⣿⣿⣿⣿⣿⣿⡀⠱⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣿⣿⣿⣿⣿⣿⣿⣦⠀⠱⣿⣿⣿⣿⣿⣿⣇⠀⢃⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢿⣿⣿⣿⣿⣿⣿⣷⡄⣹⣿⣿⣿⣿⣿⣿⣶⣾⣿⣶⣤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣴⣶⣿⣭⣍⡉⠙⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⢀⣠⣶⣿⣿⣿⣿⣿⣿⣿⣿⣷⣦⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀⠀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠉⠉⠛⠻⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡷⢂⣓⣶⣶⣶⣶⣤⣤⣄⣀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⣿⣿⠟⢀⣴⢿⣿⣿⣿⠟⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠛⠋⠉⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠤⠤⠤⠤⠙⣻⣿⣿⣿⣿⣿⣿⣾⣿⣿⡏⣠⠟⡉⣾⣿⣿⠋⡠⠊⣿⡟⣹⣿⢿⣿⣿⣿⠿⠛⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣤⣶⣤⣭⣤⣼⣿⢛⣿⣿⣿⣿⣻⣿⣿⠇⠐⢀⣿⣿⡷⠋⠀⢠⣿⣺⣿⣿⢺⣿⣋⣉⣉⣩⣴⣶⣤⣤⣄⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠛⠻⠿⣿⣿⣿⣇⢻⣿⣿⡿⠿⣿⣯⡀⠀⢸⣿⠋⢀⣠⣶⠿⠿⢿⡿⠈⣾⣿⣿⣿⣿⡿⠿⠛⠋⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠻⢧⡸⣿⣿⣿⠀⠃⠻⠟⢦⢾⢣⠶⠿⠏⠀⠰⠀⣼⡇⣸⣿⣿⠟⠉⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣴⣾⣶⣽⣿⡟⠓⠒⠀⠀⡀⠀⠠⠤⠬⠉⠁⣰⣥⣾⣿⣿⣶⣶⣷⡶⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠉⠹⠟⣿⣿⡄⠀⠀⠠⡇⠀⠀⠀⠀⠀⢠⡟⠛⠛⠋⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⠋⠹⣷⣄⠀⠐⣊⣀⠀⠀⢀⡴⠁⠣⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣤⣀⠤⠊⢁⡸⠀⣆⠹⣿⣧⣀⠀⠀⡠⠖⡑⠁⠀⠀⠀⠑⢄⣀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⣦⣶⣿⣿⣟⣁⣤⣾⠟⠁⢀⣿⣆⠹⡆⠻⣿⠉⢀⠜⡰⠀⠀⠈⠑⢦⡀⠈⢾⠑⡾⠲⣄⠀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⠀⠀⠀⣀⣤⣶⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠖⠒⠚⠛⠛⠢⠽⢄⣘⣤⡎⠠⠿⠂⠀⠠⠴⠶⢉⡭⠃⢸⠃⠀⣿⣿⣿⠡⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀
   ⠀⠀⠀⠀⠀⡤⠶⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣋⠁⠀⠀⠀⠀⠀⢹⡇⠀⠀⠀⠀⠒⠢⣤⠔⠁⠀⢀⡏⠀⠀⢸⣿⣿⠀⢻⡟⠑⠢⢄⡀⠀⠀⠀⠀
   ⠀⠀⠀⠀⢸⠀⠀⠀⡀⠉⠛⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣄⣀⣀⡀⠀⢸⣷⡀⣀⣀⡠⠔⠊⠀⠀⢀⣠⡞⠀⠀⠀⢸⣿⡿⠀⠘⠀⠀⠀⠀⠈⠑⢤⠀⠀
   ⠀⠀⢀⣴⣿⡀⠀⠀⡇⠀⠀⠀⠈⣿⣿⣿⣿⣿⣿⣿⣿⣝⡛⠿⢿⣷⣦⣄⡀⠈⠉⠉⠁⠀⠀⠀⢀⣠⣴⣾⣿⡿⠁⠀⠀⠀⢸⡿⠁⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀
   ⠀⢀⣾⣿⣿⡇⠀⢰⣷⠀⢀⠀⠀⢹⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣦⣭⣍⣉⣉⠀⢀⣀⣤⣶⣾⣿⣿⣿⢿⠿⠁⠀⠀⠀⠀⠘⠀⠀⠀⠀⠀⠀⠀⠀⠀⡰⠉⢦⠀
   ⢀⣼⣿⣿⡿⢱⠀⢸⣿⡀⢸⣧⡀⠀⢿⣿⣿⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡭⠖⠁⠀⡠⠂⠀⠀⠀⠀⠀⠀⠀⠀⢠⠀⠀⠀⢠⠃⠀⠈⣀
   ⢸⣿⣿⣿⡇⠀⢧⢸⣿⣇⢸⣿⣷⡀⠈⣿⣿⣇⠈⠛⢿⣿⣿⣿⣿⣿⣿⠿⠿⠿⠿⠿⠿⠟⡻⠟⠉⠀⠀⡠⠊⠀⢠⠀⠀⠀⠀⠀⠀⠀⠀⣾⡄⠀⢠⣿⠔⠁⠀⢸
   ⠈⣿⣿⣿⣷⡀⠀⢻⣿⣿⡜⣿⣿⣷⡀⠈⢿⣿⡄⠀⠀⠈⠛⠿⣿⣿⣿⣷⣶⣶⣶⡶⠖⠉⠀⣀⣤⡶⠋⠀⣠⣶⡏⠀⠀⠀⠀⠀⠀⠀⢰⣿⣧⣶⣿⣿⠖⡠⠖⠁
   ⠀⣿⣿⣷⣌⡛⠶⣼⣿⣿⣷⣿⣿⣿⣿⡄⠈⢻⣷⠀⣄⡀⠀⠀⠀⠈⠉⠛⠛⠛⠁⣀⣤⣶⣾⠟⠋⠀⣠⣾⣿⡟⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣿⠷⠊⠀⢰⠀
   ⢰⣿⣿⠀⠈⢉⡶⢿⣿⣿⣿⣿⣿⣿⣿⣿⣆⠀⠙⢇⠈⢿⣶⣦⣤⣀⣀⣠⣤⣶⣿⣿⡿⠛⠁⢀⣤⣾⣿⣿⡿⠁⠀⠀⠀⠀⠀⠀⠀⣸⣿⡿⠿⠋⠙⠒⠄⠀⠉⡄
   ⣿⣿⡏⠀⠀⠁⠀⠀⠀⠉⠉⠙⢻⣿⣿⣿⣿⣷⡀⠀⠀⠀⠻⣿⣿⣿⣿⣿⠿⠿⠛⠁⠀⣀⣴⣿⣿⣿⣿⠟⠀⠀⠀⠀⠀⠀⠀⠀⢠⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠰
   ====================================================
   BIOCORE@CRG Master of Pores 4. Preprocessing - N F  ~  version 4.0
   ====================================================
   
   
   Input
   ----------------------------------------------------
   pod5                      : 
   fastq                     : /Users/lcozzuto/ooo/master_of_pores/mop_preprocess/../data/fastq/*.fq.gz
   
   Reference
   ----------------------------------------------------
   reference                  : /Users/lcozzuto/ooo/master_of_pores/mop_preprocess/../anno/yeast_rRNA_ref.fa.gz
   annotation                 : 
   ref_type                   : transcriptome
   
   Output
   ----------------------------------------------------
   output                    : ./outfolder/
   email                     : 
   slackhook                 : 
   
   Actions
   ----------------------------------------------------
   basecalling               : NO
   demultiplexing            : dorado
   demulti_pod5              : ON
   filtering                 : nanoq
   mapping                   : minimap2
   counting                  : nanocount
   discovery                 : NO
   cram_conv                 : YES
   subsampling_cram          : 50
   
   Advanced
   ----------------------------------------------------
   granularity               : 1
   barcodes                  : 
   GPU                       : OFF
   
   ====================================================
   
   
   ----------------------CHECK TOOLS -----------------------------
   > basecalling will be skipped
   > demultiplexing will be skipped
   mapping : minimap2
   filtering : nanoq
   counting : nanocount
   > discovery will be skipped
   --------------------------------------------------------------
   WARN: Include with `addParams()` is deprecated -- pass params as a workflow or process input instead
   Skipping the email
   
   executor >  local (24)
   [a3/c6f206] MAPPING_MOP:ALIGN:MINIMAP2:map (wt)           [100%] 2 of 2 ✔
   [aa/f74040] SAMTOOLS_SORT:sortAln (mod)                   [100%] 2 of 2 ✔
   [6d/8d015f] SAMTOOLS_INDEX:indexBam (mod)                 [100%] 2 of 2 ✔
   [13/7306b7] bam2stats (mod)                               [100%] 2 of 2 ✔
   [22/86891e] joinAlnStats (joining aln stats)              [100%] 1 of 1 ✔
   [53/e313a4] NANOSTAT_QC:nanoStat (mod)                    [100%] 2 of 2 ✔
   [e5/b3465a] checkRef (Checking yeast_rRNA_ref.fa.gz)      [100%] 1 of 1 ✔
   [2c/328c32] bam2Cram (mod)                                [100%] 2 of 2 ✔
   [91/58b424] NANOQ_REPORT:report (mod)                     [100%] 2 of 2 ✔
   [d6/ad8c41] COUNTING:NANOCOUNT:nanoCount (mod)            [100%] 2 of 2 ✔
   [24/7aa36d] COUNTING:AssignReads (mod)                    [100%] 2 of 2 ✔
   [0f/5a98e2] COUNTING:countStats (mod)                     [100%] 2 of 2 ✔
   [9e/fc64d5] COUNTING:joinCountStats (joining count stats) [100%] 1 of 1 ✔
   [2f/23940b] MULTIQC:makeReport                            [100%] 1 of 1 ✔
   ---------------------------------------------------
               *Pipeline MOP4 completed!*             
   ---------------------------------------------------
   - Launched by `lcozzuto`
   - Started at 2025-04-10 17:07:31
   - Finished at 2025-04-10 17:07:47
   - Time elapsed: 16.2s
   - Execution status: OK
   ```nextflow run mop_preprocess.nf -with-docker -params-file params.yaml```
   ---------------------------------------------------
   
   
.. Note:: 
   The latest versions of Nextflow show a future deprecation of `addParams()`. For now just ignore this warning.
.. code-block:: console
   WARN: Include with `addParams()` is deprecated -- pass params as a workflow or process input instead

The output folders will be in `outfolder` as indicated by the parameter `output`. Inside, you have the following list of directories:

- alignment: sorted bam files and their indexes. 
- assigned: tabular file with index id and assigned chromosome or transcript
- counts: read counts per feature (transcript or gene)
- cram_files: sorted, subsampled cram files and their indexes.  
- report: multiq report

Starting from pod5 (linux or mac local)
======================







