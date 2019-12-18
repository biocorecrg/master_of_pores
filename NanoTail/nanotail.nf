#!/usr/bin/env nextflow


/* 
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '0.1'

params.help            = false
params.resume          = false

log.info """

╔╦╗┌─┐┌─┐┌┬┐┌─┐┬─┐  ┌─┐┌─┐  ╔═╗╔═╗╦═╗╔═╗╔═╗
║║║├─┤└─┐ │ ├┤ ├┬┘  │ │├┤   ╠═╝║ ║╠╦╝║╣ ╚═╗
╩ ╩┴ ┴└─┘ ┴ └─┘┴└─  └─┘└    ╩  ╚═╝╩╚═╚═╝╚═╝
                                                                                       
====================================================
BIOCORE@CRG NanoTail. Detection of polyA length (RNA) - N F  ~  version ${version}
====================================================

*****************   Input files    *********************
input_folders            : ${params.input_folders}

******* reference has to be the transcriptome **********
reference                : ${params.reference}
output                   : ${params.output}

********* nanopolish and tailfindr cmd options **********
nanopolish_opt           : ${params.nanopolish_opt} 
tailfindr_opt            : ${params.tailfindr_opt} 

email                     : ${params.email}

"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
config_report = file("$baseDir/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("$baseDir/../docs/logo_small.png")

nanopolish_opt 	   = params.nanopolish_opt
tailfindr_opt      = params.tailfindr_opt

// Output folders
outputTailfindr     = "${params.output}/Tailfindr"
outputNanoPolish    = "${params.output}/NanoPolish"
outputFinal         = "${params.output}/PolyA_final"

/*
 * Creates the channels that emits input data
 */
Channel
    .fromFilePairs( params.input_folders, type: 'dir', size: 1) 
    .ifEmpty { error "Cannot find any folder matching: ${params.input_folders}" }
    .into{folders_for_tailfindr; folders_for_nanopolish}
    
folders_for_nanopolish.map {  
        sampleID = it[0]
        folderPath = it[1][0]
        fast5 = file("${folderPath}/fast5_files")
        fastq = file("${folderPath}/fastq_files/*")
        alignment = file("${folderPath}/alignment/*")
        [ sampleID, fast5, fastq, alignment ]
    }.set{data_for_nanopolish}

folders_for_tailfindr.map {  
        sampleID = it[0]
        folderPath = it[1][0]
        fast5 = file("${folderPath}/fast5_files")
        fastq = file("${folderPath}/fastq_files/*")
        alignment = file("${folderPath}/alignment/*")
        [ sampleID, fast5, fastq ]
    }.set{data_for_tailfindr}

/*
* Estimate polyA tail size with tailfindr
*/

process tailfindr {
	publishDir outputTailfindr, pattern: "*_findr.csv",  mode: 'copy'
	tag { sampleID }  
	label 'big_mem_cpus_time'
	
	input:
	set val(sampleID), file(fast5), file(fastq) from data_for_tailfindr

	output:
	set val(sampleID), file("${sampleID}_findr.csv") into tailfindr_res

	script:
	"""
	R --slave -e "library(tailfindr); find_tails(fast5_dir = './fast5_files' , save_dir = './', ${params.tailfindr_opt}, csv_filename = \'${sampleID}_findr.csv\', num_cores = ${task.cpus})"
	"""
}

/*
* Estimate polyA tail size with nanopolish
*/
process tail_nanopolish {
	publishDir outputNanoPolish, pattern: "*.polya.estimation.tsv",  mode: 'copy'
	tag { sampleID }  
	label 'big_mem_cpus_time'
	
	input:
	file(reference)
	set val(sampleID), file(fast5), file(fastq), file(alignment) from data_for_nanopolish

	output:
	set val(sampleID), file("${sampleID}.polya.estimation.tsv") into nanopol_res

	script:
	"""
    if [ `echo ${reference} | grep ".gz"` ]; then 
        zcat ${reference} > `basename ${reference} .gz`
        myreference=`basename ${reference} .gz`
    else myreference=${reference}
    fi
    #to keep only mapped reads and remove secondary alignments 
    samtools view -bF 260 ${alignment} > ${sampleID}.bam 
	samtools index ${sampleID}.bam 
	#index reads
	nanopolish index -d ./ ${fastq}
	# polya length estimation
	nanopolish polya -r ${fastq} ${params.nanopolish_opt} -g \$myreference -t ${task.cpus} -b ${sampleID}.bam > ${sampleID}.polya.estimation.tsv
	rm \$myreference
	"""
} 

/*
process join_results {
	publishDir outputFinal,  mode: 'copy'
	tag { sampleID }  
	label 'big_mem_cpus_time'
	
	input:
	set val(sampleID), file(tail_res), file(nano_res) from tailfindr_res.join(nanopol_res)
	
	output:

	script:
	"""
	"""
	
}


/*
* Merge the results
*/



/*
* functions
*/
// Get the name from the folder
def getFolderName(sample) {
   folder_info = sample.toString().tokenize("/")
   return folder_info[-2]
}

// make named pipe 
def unzipBash(filename) { 
    cmd = filename.toString()
    if (cmd[-3..-1] == ".gz") {
    	cmd = "<(zcat ${filename})"
    }
    return cmd
}



/*
*  Finish message
*/
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

/*
 * Mail notification


workflow.onComplete {
    def subject = 'Master of Pore execution'
    def recipient = "${params.email}"
    def attachment = "${outputMultiQC}/multiqc_report.html"

    ['mail', '-s', subject, '-a', attachment, recipient].execute() << """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}
*/