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

******* reference has to be the genome **********
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
joinScript = file("$baseDir/bin/join.r")


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
    .into{folders_for_fastq; folders_for_analyisis; folders_for_genes; folders_for_bam_filtering}
    
folders_for_analyisis.map {  
        [ it[0], file("${it[1][0]}/fast5_files/*") ]
    }.transpose().into{data_for_tailfindr; data_for_nanopolish_raw}

folders_for_bam_filtering.map {  
        [ it[0], file("${it[1][0]}/alignment/*") ]
    }.set{data_for_bam_filtering}

folders_for_fastq.map {
        [ it[0], file("${it[1][0]}/fastq_files/*") ]
    }.set{fastq_data}
    
folders_for_genes.map {  
        [ it[0], file("${it[1][0]}/assigned/*") ]
    }.set{genes_for_final}


/*
* Check reference file
*/
process check_reference {
	
	input:
	file(reference)

	output:
    file("reference_sequences.fa") into checked_ref

	script:
	"""
    if [ `echo ${reference} | grep ".gz"` ]; then 
        zcat ${reference} > reference_sequences.fa
    else ln -s ${reference} reference_sequences.fa
    fi
	"""
}

/*
* Estimate polyA tail size with tailfindr
*/

process tailfindr {
	tag "${sampleID}-${fast5.simpleName}"  
	label 'big_mem_cpus'
	
	input:
	set val(sampleID), file(fast5) from data_for_tailfindr

	output:
	set val(sampleID), file("*_findr.csv") into tailfindr_res

	script:
	def fast5_index=fast5.getSimpleName()
	"""
	R --slave -e "library(tailfindr); find_tails(fast5_dir = './' , save_dir = './', ${params.tailfindr_opt}, csv_filename = \'${sampleID}-${fast5_index}_findr.csv\', num_cores = ${task.cpus})"
	"""
}

/*
* Filter bams
*/
process filter_bam {
	tag { sampleID }  
	label 'big_mem_cpus'
	
	input:
	file(reference)
	set val(sampleID), file(alignment) from data_for_bam_filtering

	output:
	set val(sampleID), file("${sampleID}_filt.bam"), file("${sampleID}_filt.bam.bai") into filt_bam_for_nanopolish, filt_bam_for_intersect

	script:
	"""
    #to keep only mapped reads and remove secondary alignments 
    samtools view -@ {task.cpus} -bF 260 ${alignment} > ${sampleID}_filt.bam 
	samtools index ${sampleID}_filt.bam 
	"""
}


data_for_nanopolish_raw.combine(filt_bam_for_nanopolish, by:0).combine(fastq_data, by:0)
.set{data_for_nanopolish}



/*
* Estimate polyA tail size with nanopolish
*/
process tail_nanopolish {
	tag "${sampleID}-${fast5.simpleName}"  
	label 'big_mem_cpus_time'
	
	input:
	file(checked_ref)
	set val(sampleID), file(fast5), file(alignment), file(alnindex), file(fastq) from data_for_nanopolish

	output:
	set val(sampleID), file("*.polya.estimation.tsv") into nanopol_res

	script:
	def fast5_index=fast5.getSimpleName()
	"""
	#index reads
	nanopolish index -d ./ ${fastq}
	# polya length estimation
	nanopolish polya -r ${fastq} ${params.nanopolish_opt} -g ${checked_ref} -t ${task.cpus} -b ${alignment} > ${sampleID}-${fast5_index}.polya.estimation.tsv
	"""
} 

process collect_nanopolish_results {
	publishDir outputNanoPolish, pattern: "*.polya.estimation.tsv",  mode: 'copy'
	tag { sampleID }  
	
	input:
	set val(sampleID), file("nanopol_*") from nanopol_res.groupTuple()
	
	output:
    set val(sampleID), file("${sampleID}.nanopol.len") into nanopol_len
    file ("*.polya.estimation.tsv") 

	script:
	"""
	head -n 1 nanopol_1 >> ${sampleID}.polya.estimation.tsv 
	grep --no-filename  -v "leader_start" nanopol_* | grep -v "READ_FAILED_LOAD" >> ${sampleID}.polya.estimation.tsv 
	awk -F"\t" '{if (\$10=="PASS") print \$1"\t"\$9}' ${sampleID}.polya.estimation.tsv > ${sampleID}.nanopol.len
	"""

}

process collect_tailfindr_results {
	publishDir outputTailfindr, pattern: "*_findr.csv",  mode: 'copy'
	tag { sampleID }  
	
	input:
	set val(sampleID), file("tailfin_*") from tailfindr_res.groupTuple()
	
	output:
    set val(sampleID), file("${sampleID}.findr.len") into tailfindr_len
    file ("*_findr.csv") 

	script:
	"""
	head -n 1 tailfin_1 >> ${sampleID}_findr.csv
	grep --no-filename  -v "tail_start" tailfin_* >> ${sampleID}_findr.csv
	awk -F"," '{if (\$5!="" && \$1!="read_id") print \$1"\t"\$5}' ${sampleID}_findr.csv > ${sampleID}.findr.len
	"""


}




/*
*/
process join_results {
	publishDir outputFinal,  mode: 'copy'
	tag { sampleID }  
	
	input:
	set val(sampleID), file(nanopol), file(tailfindr), file(genes) from nanopol_len.join(tailfindr_len).join(genes_for_final)
	file(joinScript)
	
	output:
	file("${sampleID}_*")
	
	script:
	"""
	Rscript --vanilla join.r ${tailfindr} ${nanopol} ${genes} ${sampleID}
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
    println "Pipeline BIOCORE@CRG Master of Pore completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
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
