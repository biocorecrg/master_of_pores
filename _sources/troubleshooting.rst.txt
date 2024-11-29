.. _home-page-troubleshooting:

*****************
Troubleshooting
*****************

.. autosummary::
   :toctree: generated


Failure in downloading singularity images
================================================================
Sometimes we have problems in downloading the singularity images. You can pre-download them by using the script **get_singularity.py** in **MOP4/BioNextflow/scripts**. You need to inspect your pipeline and store the results as here:


.. code-block:: console

   nextflow inspect -params-file params.yaml -profile slurm mop_preprocess.nf > inspect.json

Then you can feed the json file to **get_singularity.py** and indicate where to store the images:

.. code-block:: console

   python get_singularity.py -j inspect.json -c ./singularity_folder


Memory failures in mop_preprocess and mop_mod
================================================================

Sometimes FastQC, Epinano, or other tools can run out of memory because the input data contain low-quality reads. It is always a good idea to double-check it and eventually filter them per size / quality using either **nanoq** or **nanofilt**.
