// Output folders

//mop_preprocess
outputFastq    = "${params.output}/fastq_files"
outputFast5    = "${params.output}/pod5_files"
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
outputF5CFlow        = "${params.output}/f5c_flow"
outputM6Anet         = "${params.output}/m6anet_flow"
outputModPhredFlow   = "${params.output}/modphred_flow"
outputModKitFlow     = "${params.output}/modkit_flow"

// mop_dna
outputClairS         = "${params.output}/clair_flow"
outputSniffles       = "${params.output}/sniffles_flow"
outputAnnotation     = "${params.output}/annotation"
outputFiltered       = "${params.output}/filtered_vcf"
