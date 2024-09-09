#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* 
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '2.0'

params.help            = false
params.resume          = false

log.info """

╔╦╗╔═╗╔═╗  ╔╦╗╔═╗╦╦  
║║║║ ║╠═╝   ║ ╠═╣║║  
╩ ╩╚═╝╩     ╩ ╩ ╩╩╩═╝
                                                                                       
====================================================
BIOCORE@CRG Master of Pores 2. Estimating PolyA tail - N F  ~  version ${version}
====================================================

*****************   Input files    *********************
input_path               : ${params.input_path}
output                   : ${params.output}
pars_tools				 : ${params.pars_tools}

******* reference has to be the genome **********
reference                : ${params.reference}

email                     : ${params.email}

************************* Flows *******************************
tailfindr                             	: ${params.tailfindr}
nanopolish                              : ${params.nanopolish}

email                                   : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// include functions, outdirs from other files
evaluate(new File("../outdirs.nf"))
def local_modules = file("$baseDir/../local_modules.nf")
def subworkflowsDir = "${baseDir}/../BioNextflow/subworkflows"
joinScript = file("$baseDir/bin/join.r")

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"

def flows = [:]
flows["tailfindr"] = params.tailfindr
flows["nanopolish"] = params.nanopolish

include { getParameters; checkRef } from "${local_modules}" 

progPars = getParameters(params.pars_tools)

include { ESTIMATE_TAIL as TAILFINDR_ESTIMATE_TAIL } from "${subworkflowsDir}/chem_modification/tailfindr" addParams(LABEL: 'big_cpus_retry', EXTRAPARS: progPars["tailfindr--tailfindr"])
include { GET_VERSION as TAILFINDR_VER } from "${subworkflowsDir}/chem_modification/tailfindr" 
include { GET_VERSION as SAMTOOLS_VER; INDEX as SAMTOOLS_INDEX } from "${subworkflowsDir}/misc/samtools"
include { POLYA_LEN as NANOPOLISH_POLYA_LEN } from "${subworkflowsDir}/chem_modification/nanopolish" addParams(LABEL: 'big_cpus',  OUTPUT: outputNanopolish, EXTRAPARS: progPars["nanopolish--nanopolish"])
include { GET_VERSION as NANOPOLISH_VER } from "${subworkflowsDir}/chem_modification/nanopolish" 

include { reshapeSamples } from "${local_modules}" 
include { collect_tailfindr_results} addParams(OUTPUT: outputTailFindr) from "${local_modules}"
include { join_nanotail_results } addParams(OUTPUT: outputFinalPolyA) from "${local_modules}"
include { filter_bam} addParams(LABEL: 'big_cpus') from "${local_modules}"
 

Channel.fromFilePairs("${params.input_path}/alignment/*_s.bam", size: 1).set{bams}
Channel.fromFilePairs("${params.input_path}/alignment/*_s.bam.bai", size: 1).set{bais}
Channel.fromFilePairs("${params.input_path}/fastq_files/*.fq.gz", size: 1).set{fastqs}
Channel.fromFilePairs("${params.input_path}/fast5_files/*", size: 1, type: 'dir').set{fast5_folders}
Channel.fromFilePairs("${params.input_path}/assigned/*.assigned", size: 1).set{assigned}



fast5_folders.map{
	[ it[0], file("${it[1][0]}/*.fast5") ]
}.transpose().set{fast5_files_4_np}

fast5_files_4_np.map{
	[ "${it[0]}---${it[1].simpleName}", it[1] ]
}.set{fast5_files_4_tf}


workflow {	

	if (params.tailfindr == "YES") {
		tail_estim = TAILFINDR_ESTIMATE_TAIL(fast5_files_4_np)
		//resh_tail_estim = reshapeSamples(tail_estim)
		//tail_estim.groupTuple().view()
		tailres = collect_tailfindr_results(tail_estim.groupTuple())
	}
	if (params.nanopolish == "YES") {
		ref_file = checkRef(reference)
		filt_bams = filter_bam(ref_file, bams)
		filt_bais = SAMTOOLS_INDEX(filt_bams)
		nanores = NANOPOLISH_POLYA_LEN(fast5_files_4_np, bams, bais, fastqs, ref_file) 
	}
	if (params.tailfindr == "YES" && params.nanopolish == "YES") {
                log.info "Joining results"
				//nanores.filtered_est.view()
                join_nanotail_results(nanores.filtered_est.join(tailres.length).join(assigned), joinScript)

	}

	all_ver = TAILFINDR_VER().mix(NANOPOLISH_VER())
	.mix(SAMTOOLS_VER())
	.collectFile(name: 'tool_version.txt', newLine: false, storeDir:outputFinalPolyA)

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
*/

if (params.email == "yourmail@yourdomain" || params.email == "") { 
    log.info 'Skipping the email\n'
}
else {
    log.info "Sending the email to ${params.email}\n"

    workflow.onComplete {

    def msg = """\
        Pipeline BIOCORE@CRG Master of Pore 2 modification module's execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()

        sendMail(to: params.email, subject: "Master of Pore 2 execution", body: msg)
    }
}
