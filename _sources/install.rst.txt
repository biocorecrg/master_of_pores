.. _home-page-install:

**************
Get Started
**************

.. autosummary::
   :toctree: generated

Please install nextflow `Nextflow <https://www.nextflow.io/>`_ and either `Singularity <https://sylabs.io/>`_ or `Docker <https://www.docker.com/>`_ before.

For installing Nextflow you need a POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and Java 11 (or later, up to 17). Windows system is supported through WSL. For the installation of Nextflow just run:

.. code-block:: console

  curl -s https://get.nextflow.io | bash

To install the pipeline you need to download the repo:

.. code-block:: console

   git clone --depth 1 --recurse-submodules https://github.com/biocorecrg/MOP4.git


Testing
============

.. code-block:: console

  cd mop_preprocess

  nextflow run mop_preprocess.nf -params-file params.f5.yaml -with-singularity -bg -profile local > log

.. tip::

  You can replace ```-with-singularity``` with ```-with-docker``` if you want to use the docker engine.

Profiles
============
Some nextflow configuration files are stored within the folder **conf** and can be selected using different profiles. Currently, we have:

- ci:              for continuous integration testing (low resources)
- local:           for being used in a laptop without GPU support
- m1mac:           for running the containers in emulation for being used on M1/M2/M3 Apple processors.
- sge:             for being used in an HPC with Sun Grid Engine
- cluster or crg:  for being used in the custom HPC environment at CRG
- slurm:           for being used in an HPC with SLURM
- awsbatch:        for being used in Amazon AWS cloud infrastructure
- clean:           for removing intermediate files (resume won't work)

Experimental cleanup
======================
The nextflow pipelines can occupy a lot of space in the work folders. If you are not interested in caching and resuming you might want to use the clean profile in addition to the one you choose for the infrastructure.

.. code-block:: console

  cd mop_preprocess

  nextflow run mop_preprocess.nf -params-file params.f5.yaml -with-singularity -bg -profile local,clean > log

Specify the place for singularity image download
===============
Setting the following variables in your .bashrc or .bash_profile will specify the place for downloading the images:

.. code-block:: console
  export APPTAINERENV_TMPDIR="MYPATH"
  export APPTAINERENV_NXF_TASK_WORKDIR="MYPATH"
  export NXF_APPTAINERENV_LIBRARYDIR="MYPATH"
  export NXF_SINGULARITY_CACHEDIR="MYPATH"

Running on slurm HPC by submitting the nextflow job 
===============
You can use the script `launch_nf.sh`for submitting the nextflow jobs as follows:


.. code-block:: console

   sbatch launch_nf.sh nextflow run run mop_preprocess.nf -params-file params.f5.yaml -ansi-log false -with-singularity -profile newcrg 



