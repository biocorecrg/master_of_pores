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


Error with demultiplexing with Dorado
================================================================

When running dorado for demultiplexing you would need to add the option `--no-trim` or you could get errors. This is because of this `issue <https://github.com/nanoporetech/dorado/issues/539>`_:

