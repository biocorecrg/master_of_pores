.. _home-page-install:

**************
How to install
**************

.. autosummary::
   :toctree: generated

Please install nextflow `Nextflow <https://www.nextflow.io/>`_ and either `Singularity <https://sylabs.io/>`_ or `Docker <https://www.docker.com/>`_ before.

For installing Nextflow you need a POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and Java 8 (or later, up to 17). Windows system is supported through WSL. For the installation of Nextflow just run:

.. code-block:: console

  curl -s https://get.nextflow.io | bash

For installing the pipeline you need to download the repo:

.. code-block:: console

  git clone --depth 1 --recurse-submodules https://github.com/biocorecrg/MOP2.git


You can use **INSTALL.sh** to download the **guppy 3.4.5** or you can download the version you prefer by adding an extra parameter. 

.. note::
  
  Please consider that the support of VBZ compression of fast5 started with version 3.4.X. 


.. code-block:: console
  
  cd MoP2; bash INSTALL.sh 

or

.. code-block:: console

  cd MoP2; bash INSTALL.sh 4.0.15
 
 
Testing
============

.. code-block:: console

  cd mop_preprocess

  nextflow run mop_preprocess.nf -with-singularity -bg -profile local > log

.. tip::

  You can replace ```-with-singularity``` with ```-with-docker``` if you want to use the docker engine.


Apple M1 processor
====================

Use the profile **m1mac** for running on machines with Apple M1 processor.


.. code-block:: console

  cd mop_preprocess

  nextflow run mop_preprocess.nf -with-singularity -bg -profile m1mac > log


