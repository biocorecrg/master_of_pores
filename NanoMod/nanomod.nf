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
BIOCORE@CRG NanoDirectRNA. Detection of modification and polyA length (RNA) - N F  ~  version ${version}
====================================================

*****************   Input files    *******************
input_path                              : ${params.input_path}
comparison                              : ${params.comparison}

********** reference has to be the genome *************
reference                               : ${params.reference}
output                                  : ${params.output}

coverage                                : ${params.coverage}
************* tombo and epinano params ****************
tombo_opt                               : ${params.tombo_opt}

epinano_opt                             : ${params.epinano_opt}

email                                   : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
//config_report = file("$baseDir/config.yaml")
//if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("$baseDir/../docs/logo_small.png")

model_folder = file("$baseDir/models/")
if( !model_folder.exists() ) exit 1, "Missing folders with EpiNano's models!"
joinScript = file("$baseDir/bin/join.r")

tombo_opt    	= params.tombo_opt
epinano_opt     = params.epinano_opt

// Output folders
outputtombo   = "${params.output}/Tombo"
outputEpinano      = "${params.output}/Epinano"
//outputReport   = file("${outputMultiQC}/multiqc_report.html")

/*
* move old multiQCreport

if( outputReport.exists() ) {
  log.info "Moving old report to multiqc_report.html multiqc_report.html.old"
  outputReport.moveTo("${outputMultiQC}/multiqc_report.html.old")
}
*/

compfile = file(params.comparison)
if( !compfile.exists() ) exit 1, "Missing comparison file: ${compfile}. Specify path with --comparisons"

/*
 * Creates the channels with comparisons
 */
 Channel
    .from(compfile.readLines())
    .map { line ->
        list = line.split("\t")
        if (list[0]!= "") {
            def sampleID = list[0]
            def ctrlID = list[1]
            [ sampleID, ctrlID ]
        }
    }
    .into{ id_to_tombo_fast5; id_to_tombo_idx; luca; id_to_epinano; id_for_resquiggling}


/*
 * Creates the channels with samples and KOs for EpiNano
 */
 Channel
    .from(compfile.readLines())
    .map { line ->
        list = line.split("\t")
        if (list[0]!= "") {
            list[0] 
        }
    }
    .set{ samples_for_epinano_filtering }

 Channel
    .from(compfile.readLines())
    .map { line ->
        list = line.split("\t")
        if (list[1]!= "") {
            list[1] 
        }
    }
    .set{ ko_for_epinano_filtering }

id_for_resquiggling.flatten().unique().map {
    [it, file("${params.input_path}/${it}/fast5_files/*.fast5")]
}.transpose().set{fast5_for_resquiggling}

id_to_epinano.flatten().unique().map {
    [it, file("${params.input_path}/${it}/alignment/*.bam")]
}.transpose().set{bams_for_variant_calling}



/*
* Perform resquiggling
*/

process resquiggling {
    tag {"${idsample}-${fast5.simpleName}"}  
    label 'big_mem_cpus'
	
    input:
    set idsample, file(fast5) from fast5_for_resquiggling
    file(reference)
    
    output:
    file ("*.resquiggle.failed.tsv") into failed_resquiggles
    set idsample, file ("${idsample}-${fast5.simpleName}") into singlefast5
    set idsample, file (".*.tombo.index") into tombo_indexes

    script:
    def infolder = "${idsample}-${fast5.simpleName}"
    
    """ 
	#from multifast5 to singlefast5
    mkdir ${infolder};
    multi_to_single_fast5 -i ./ -s ./ -t ${task.cpus}; 
    rm ./filename_mapping.txt; 
    mv ./*/*.fast5 ${infolder};
    # resquiggling
    tombo resquiggle ${infolder} ${reference} --rna --processes ${task.cpus} --overwrite --failed-reads-filename ${infolder}.resquiggle.failed.tsv 
    """
}

/*
* Group data together with indexes
*/
singlefast5.groupTuple().into{grouped_single_fast5_A; grouped_single_fast5_B}
id_to_tombo_fast5.combine(grouped_single_fast5_A, by: 0).map {
	[ it[1], it[0], it[2] ]
}.combine(grouped_single_fast5_B, by: 0)map {
	[ "${it[1]}--${it[0]}", it[2], it[3]]
}.set{fast5_for_tombo_modifications}

tombo_indexes.groupTuple().into{grouped_indexes_A; grouped_indexes_B}
id_to_tombo_idx.combine(grouped_indexes_A, by: 0).map {
	[ it[1], it[0], it[2] ]
}.combine(grouped_indexes_B, by: 0)map {
	[ "${it[1]}--${it[0]}", it[2], it[3] ]
}.set{idx_for_tombo_modifications}

fast5_for_tombo_modifications.combine(idx_for_tombo_modifications, by: 0).set{data_for_tombo_modifications}

/*
* detect modification 
*/

process getModifications {
    label 'big_mem_cpus'
    tag {"${combID}"}  
	publishDir outputtombo, pattern: "*.significant_regions.fasta",  mode: 'copy'
        
    input:
    file(reference)
    set val(combID), file(fast5s_A), file(fast5s_B), file(index_A), file(index_B) from data_for_tombo_modifications
    
    output:
    file ("*.significant_regions.fasta") into sign_tombo_regions

    script:
	def reference_cmd = unzipFile(reference, "reference.fa")
	def folder_names = "${combID}".split("--")
	def folder_name_A = folder_names[0]
	def folder_name_B = folder_names[1]
	"""
	${reference_cmd}
	mkdir ${folder_name_A} ${folder_name_B}
	mv ${fast5s_A} ${folder_name_A}
	mv ${index_A} ${folder_name_A}
	mv ${fast5s_B} ${folder_name_B}
	mv ${index_B} ${folder_name_B}
	tombo detect_modifications model_sample_compare --minimum-test-reads ${params.coverage} --fast5-basedirs ${folder_name_A}/* --control-fast5-basedirs ${folder_name_B}/* --statistics-file-basename ${folder_name_A}_${folder_name_B}_model_sample_compare --rna --per-read-statistics-basename ${folder_name_A}_${folder_name_B}_per-read-statistics --processes ${task.cpus}
	tombo text_output signif_sequence_context ${tombo_opt} --num-regions 1000000000 --statistics-filename ${folder_name_A}_${folder_name_B}_model_sample_compare.tombo.stats  --genome-fasta reference.fa --fast5-basedirs ${folder_name_A} --sequences-filename ${folder_name_A}_${folder_name_B}.significant_regions.fasta 
	rm reference.fa
	"""
}


/*
* Perform preprocessing for Epinano
*/

process index_reference {

    input:
    file(reference)
    
    output:
    set file("reference.fa"), file("*.dict"), file ("*.fai") into indexes
    
    script:
	"""
	if [ `echo ${reference} | grep ".gz"` ]; then 
   		zcat ${reference} > reference.fa
	else 
        ln -s ${reference} reference.fa
	fi
	\$PICARD CreateSequenceDictionary R=reference.fa O=reference.dict
	samtools faidx reference.fa
	"""
}


/*
* Calling variants for Epinano
*/

process CallVariantsForEpinano {

    tag {"${sampleID}"}  
	
    input:
    set val(sampleID), file(alnfile) from bams_for_variant_calling
    set file(reference), file(dict_index), file(faiidx) from indexes

    output:
    set sampleID, file("${sampleID}.tsv") into variants_for_frequency
   
    script:
	"""
	samtools view -h ${alnfile} -F 256 | \$SAM2TSV -r ${reference} >${sampleID}.tsv
	"""
}

/*
* Calc variant frequencies for Epinano
*/

process calcVarFrequenciesForEpinano {
	publishDir outputEpinano,  pattern: "*.csv.gz", mode: 'copy'

    tag {"${sampleID}"}  
    label 'big_mem_cpus'
	
    input:
    set val(sampleID), file(tsvfile) from variants_for_frequency
    
    output:
    set val(sampleID), file("*per_site_var.5mer.csv.gz") into per_site_vars
    file("*.csv.gz")
       
    script:
	"""
	TSV_to_Variants_Freq.py3 -f ${tsvfile} -t ${task.cpus}
	for i in *.csv; do gzip \$i; done
	"""
}


//per_site_vars.map{
//	[ it[0], it[1].splitText( by: 1, keepHeader:true, file:true ) ]
//}.transpose().set{per_site_vars_splitted}















	








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

def unzipFile(zipfile, filename) {
    cmd = "cp ${zipfile} ${filename}"
    filestring = zipfile.toString()
    if (filestring[-3..-1] == ".gz") {
    	cmd = "zcat ${zipfile} > ${filename}"
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
*/

if (params.email == "yourmail@yourdomain" || params.email == "") { 
    log.info 'Skipping the email\n'
}
else {
    log.info "Sending the email to ${params.email}\n"

    workflow.onComplete {

    def msg = """\
        NanoMod module's execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()

        sendMail(to: params.email, subject: "Master of Pore execution", body: msg)
    }
}
