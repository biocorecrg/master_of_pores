*******************
Benchmark
*******************

We tested MoP on two minION runs using the CRG's HPC where we can run up to 100 jobs in parallel (maximum 8 CPUs each) and using up to 10 GPU cards (GeForce RTX 2080 Ti). The test dataset was published at `ENA <https://www.ebi.ac.uk/>`_ with the accession `ERR5296640 <https://www.ebi.ac.uk/ena/browser/view/ERR5296640>`__  for pU samples and `ERR5303454 <https://www.ebi.ac.uk/ena/browser/view/ERR5303454>`__ for Nm samples.
 


.. list-table:: Dataset

 * - 
   - MOP_PREPROCESS
   - MOP_MOD
   - MOP_TAIL
   - MOP_CONSENSUS
 * - Input data
   - 95 Gb 
   - 137 Gb 
   - 137 Gb 
   - 14 Mb
 * - Execution time
   - 10 hours
   - 6 hours
   - 2.5 hours 
   - 3 mins
 * - Work folder
   - 382 Gb
   - 595 Gb
   - 3 Gb
   - 25 Mb
 * - Output folder
   - 137 Gb
   - 14 Mb
   - 76 Mb
   - 13 Mb
