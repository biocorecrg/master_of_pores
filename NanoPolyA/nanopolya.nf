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

fast5_sample             : ${params.fast5_sample}
fast5_control            : ${params.fast5_control}

*********** reference has to be the transcriptome ***************
reference                : ${params.reference}
output                   : ${params.output}

******* nanopolish and tailfindr are currently supported *******
tailfinder                : ${params.tailfinder} 
tailfinder_opt            : ${params.tailfinder_opt} 

email                     : ${params.email}

**************** Only required for nanopolish ******************
fastq_sample             : ${params.fastq_sample}
fastq_control            : ${params.fastq_control}


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

tailfinder 		    = params.tailfinder
tailfinder_opt      = params.tailfinder_opt

// Output folders
outputPolyA    = "${params.output}/PolyA_tail"
//outputReport   = file("${outputMultiQC}/multiqc_report.html")

/*
* move old multiQCreport

if( outputReport.exists() ) {
  log.info "Moving old report to multiqc_report.html multiqc_report.html.old"
  outputReport.moveTo("${outputMultiQC}/multiqc_report.html.old")
}
*/

/*
 * Creates the channels that emits fast5 files
 */
Channel
    .fromPath( params.fast5_sample)                                             
    .ifEmpty { error "Cannot find any sample file matching: ${params.fast5_sample}" }
    .set { fast5_sample_4_to_tail}

Channel
    .fromPath( params.fast5_control)                                             
    .ifEmpty { error "Cannot find any control file matching: ${params.fast5_control}" }
    .set { fast5_control_4_tail}

// Get the name from the folder
folder_sample_name = getFolderName(params.fast5_sample)
folder_ctrl_name = getFolderName(params.fast5_control)

// Get the channel for tail finder
data_sample_4_tails = tailFindrChannel(fast5_sample_4_to_tail, tailfinder, params.bam_sample, params.fastq_sample)
data_ctrl_4_tails = tailFindrChannel(fast5_control_4_tail, tailfinder, params.bam_control, params.fastq_control)

data_sample_4_tails.mix(data_ctrl_4_tails).into{
	data_4_tailfindr; data_4_nanopolish;
}

/*
* Estimate polyA tail size
*/

if (tailfinder == "tailfindr") {
	process tailfindr {
		publishDir outputPolyA, pattern: "*_findr.csv",  mode: 'copy'
	    tag { folder_name }  
	    label 'big_mem_cpus'
        
	    input:
	    file(reference)
	    set folder_name, file(fast5_files) from data_4_tailfindr
    
	    output:
    	file("${folder_name}_findr.csv")
    
    	script:
	    """
		R --slave -e "library(tailfindr); find_tails(fast5_dir = './' , save_dir = 'output', ${params.tailfinder}, csv_filename = \'${folder_name}_findr.csv\', num_cores = ${task.cpus})"
    	"""
	}
} else if (tailfinder == "nanopolish") {
	process tail_nanopolish {
		publishDir outputPolyA, pattern: "*.polya.estimation.tsv",  mode: 'copy'
	    tag { folder_name }  
	    label 'big_mem_cpus'
        
	    input:
	    file(reference)
	    set folder_name, file(bam_file), file(fastq_file), file(fast5_files) from data_4_nanopolish
    
	    output:
    	file("${folder_name}.polya.estimation.tsv")
    
    	script:
    	def reference_cmd = unzipBash(reference)
		"""
		#index bam
		samtools index ${bam_file}
		#index reads
		nanopolish index -d ./ ${fastq_file}
		# polya length estimation
		nanopolish polya -r ${fastq_file} ${params.tailfinder} -g ${reference_cmd} -t ${task.cpus} -b ${bam_file} > ${folder_name}.polya.estimation.tsv
		"""
	} 
} else {
     println("skipping polyA analysis")
}

/*
* functions
*/
// Get the name from the folder
def getFolderName(sample) {
   folder_info = sample.toString().tokenize("/")
   return folder_info[-2]
}

def tailFindrChannel(fast5s, tailfinder, bamfile, fastqfile) {
    fast5_4_tailfindr =Channel.empty()
    if (tailfinder == "nanopolish") {
	    fast5s.map { 
			[getFolderName(it), it]
		}.groupTuple().map { 
		[it[0], file(bamfile), file(fastqfile), it[1]]
		}.set{ fast5_4_tailfindr}
	} else {
	    fast5s.map { 
			[getFolderName(it), it]
		}.groupTuple().set{ fast5_4_tailfindr}		
	}
	return fast5_4_tailfindr
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