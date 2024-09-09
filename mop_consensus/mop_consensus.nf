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

╔╦╗╔═╗╔═╗  ╔═╗╔═╗╔╗╔╔═╗╔═╗╔╗╔╔═╗╦ ╦╔═╗
║║║║ ║╠═╝  ║  ║ ║║║║╚═╗║╣ ║║║╚═╗║ ║╚═╗
╩ ╩╚═╝╩    ╚═╝╚═╝╝╚╝╚═╝╚═╝╝╚╝╚═╝╚═╝╚═╝
                                                                                       
====================================================
BIOCORE@CRG Master of Pores 2. Get consensus modifications - N F  ~  version ${version}
====================================================

*****************   Input     *********************
input_path               : ${params.input_path}
output                   : ${params.output}
comparison               : ${params.comparison}
padsize                  : ${params.padsize}

******* reference has to be the genome **********
reference                : ${params.reference}
email                    : ${params.email}

"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// include functions, outdirs from other files
//evaluate(new File("../outdirs.nf"))
def local_modules = file("$baseDir/../local_modules.nf")
def subworkflowsDir = "${baseDir}/../BioNextflow/subworkflows"
def nanoConScript = file("$baseDir/bin/NanoConsensus.R")
def nanoScript = file("$baseDir/bin/scripts")

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"


include { checkRef; mapIDPairs; indexFasta } from "${local_modules}" 

include { nanoConsensus } from "${local_modules}" addParams(OUTPUT: params.output)

/*
 * Creates the channels with comparisons
 */

compfile = file(params.comparison)
if( !compfile.exists() ) exit 1, "Missing comparison file: ${compfile}. Specify path with --comparisons"

Channel
    .from(compfile.readLines())
    .map { line ->
        list = line.split("\t")
        if (list[0]!= "") {
            def sampleID = list[0]
            def ctrlID = list[1]
            [ sampleID, ctrlID ]
        }
    }.set {comparisons}

def padsize = (params.padsize == 0 ? 1 : params.padsize)


workflow {	
	comparisons.flatten().unique().set{unique_samples}

	unique_samples.map {
 	  	 [it, file("${params.input_path}/epinano_flow/${it}.plus_strand.per.site.csv.gz")]
	}.transpose().set{epinano}

	unique_samples.map {
 	  	 [it, file("${params.input_path}/nanopolish-compore_flow/${it}.csv.gz")]
	}.transpose().set{nanopolish}
	
	comparisons.map{
		["${it[0]}---${it[1]}", file("${params.input_path}/tombo_flow/${it[0]}---${it[1]}_lsc.plus_Tombo_Output.tsv.gz")]
	}.transpose().set{tombo}

	comparisons.map{
		["${it[0]}---${it[1]}", file("${params.input_path}/nanopolish-compore_flow/**/${it[0]}_vs_${it[1]}nanocompore_results.tsv.gz")]
	}.transpose().set{nanocomp}

	epinano_combs = mapIDPairs(comparisons, epinano).map{
		["${it[0]}---${it[1]}", it[2], it[3]]
	}
	nanopolish_combs = mapIDPairs(comparisons, nanopolish).map{
		["${it[0]}---${it[1]}", it[2], it[3]]
	}

	ref_file = checkRef(reference)
	ref_ind = indexFasta(ref_file)
	ref_ind.splitText().map{
		def vals = it.replaceAll("[\n\r]", "").split("\t")
		def chrName = vals[0]
		def chrStart = padsize
		def chrEnd = (vals[1].toInteger()-padsize)
		if (chrEnd-chrStart>(1)) {
	        [ vals[0], chrStart,  chrEnd]
	    }
	}.set{transcript_coords}
		
	data_to_process = epinano_combs.join(nanopolish_combs).join(tombo).join(nanocomp)
	//data_to_process.view()
	nanoConsensus(nanoConScript, nanoScript, ref_file, data_to_process.combine(transcript_coords))
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
