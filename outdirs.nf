// Output folders
//mop_preprocess

outputFastq    = "${params.output}/fastq_files"
outputFast5    = "${params.output}/fast5_files"
outputQual     = "${params.output}/QC_files"
outputMultiQC  = "${params.output}/report"
outputMapping  = "${params.output}/alignment"
outputCRAM     = "${params.output}/cram_files"
outputCounts   = "${params.output}/counts"
outputAssembly = "${params.output}/assembly"

//outputVars     = "${params.output}/variants"
outputAssigned = "${params.output}/assigned"

// mop_mod
outputEpinanoFlow    = "${params.output}/epinano_flow"
outputNanoPolComFlow = "${params.output}/nanopolish-compore_flow"
outputTomboFlow      = "${params.output}/tombo_flow"

// MOP_TAIL
outputTailFindr     = "${params.output}/tailfindr_flow"
outputNanopolish    = "${params.output}/nanopolish_flow"
outputFinalPolyA    = "${params.output}/polya_common"




