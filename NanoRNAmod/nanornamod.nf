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

fast5_sample             : ${params.fast5_sample}
fast5_control            : ${params.fast5_control}

*********** reference has to be the transcriptome ***************
reference                : ${params.reference}

multi5_sample            : ${params.multi5_sample}
multi5_control           : ${params.multi5_control}

output                   : ${params.output}
granularity              : ${params.granularity}

********** tombo and epinano are currently supported ************
modfinder                : ${params.modfinder}
modfinder_opt            : ${params.modfinder_opt}

******* nanopolish and tailfindr are currently supported *******
tailfinder                : ${params.tailfinder} 
tailfinder_opt            : ${params.tailfinder_opt} 

email                     : ${params.email}

*********** Only required for nanopolish / EPINANO *************
bam_sample               : ${params.bam_sample}
bam_control              : ${params.bam_control}

**************** Only required for nanopolish ******************
fastq_sample             : ${params.fastq_sample}
fastq_control            : ${params.fastq_control}


"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// check multi5. 
if (params.multi5_sample != "YES" && params.multi5_sample != "NO") exit 1, "Please specify YES or NOT in case the data is in multi fast5 format or sigle fast5"
if (params.multi5_control != "YES" && params.multi5_control != "NO") exit 1, "Please specify YES or NOT in case the data is in multi fast5 format or sigle fast5"

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
config_report = file("$baseDir/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("$baseDir/../docs/logo_small.png")

modfinder   		= params.modfinder
modfinder_opt  	    = params.modfinder_opt
tailfinder 		    = params.tailfinder
tailfinder_opt      = params.tailfinder_opt

// Output folders
outputModifs   = "${params.output}/RNA_modifs"
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
    .into { fast5_sample_4_to_batches; fast5_sample_4_to_tail; fast5_sample_4_name; fast5_sample_4_epinano}

Channel
    .fromPath( params.fast5_control)                                             
    .ifEmpty { error "Cannot find any control file matching: ${params.fast5_control}" }
    .into { fast5_control_4_to_batches; fast5_control_4_tail; fast5_control_4_name; fast5_control_4_epinano}

// Get the name from the folder
folder_sample_name = getFolderName(params.fast5_sample)
folder_ctrl_name = getFolderName(params.fast5_control)

//Make channels for modfinder
fast5_sample_4_modfinder = getModfinderChannel(fast5_sample_4_to_batches, params.granularity, folder_sample_name, params.multi5_sample)
fast5_control_4_modfinder = getModfinderChannel(fast5_control_4_to_batches, params.granularity, folder_ctrl_name, params.multi5_control)

// Get the channel for tail finder
data_sample_4_tails = tailFindrChannel(fast5_sample_4_to_tail, tailfinder, params.bam_sample, params.fastq_sample)
data_ctrl_4_tails = tailFindrChannel(fast5_control_4_tail, tailfinder, params.bam_control, params.fastq_control)

data_sample_4_tails.mix(data_ctrl_4_tails).into{
	data_4_tailfindr; data_4_nanopolish;
}

aln_4_epinano = Channel.empty()
fast5_for_epinano = Channel.empty()
if (modfinder == "epinano") {
	aln_4_epinano = Channel.from( [folder_sample_name, file(params.bam_sample)] , [folder_ctrl_name, file(params.bam_control)] )
    
    fast5_sample_4_epinano.map {
    	[folder_sample_name, file(it)]
    }.mix(fast5_control_4_epinano.map {
    	[folder_ctrl_name, file(it)]
    }).set{fast5_for_epinano}
}

/*
* Perform resquiggling
*/
process resquiggling {
    tag {"${modfinder}-${folder_name}-${idfile}"}  
    label 'big_mem_cpus'

	when:
	modfinder == "tombo"
	
    input:
    set folder_name, idfile, multi5, file(fast5) from fast5_sample_4_modfinder.mix(fast5_control_4_modfinder)
    file(reference)
    
    output:
    file ("${folder_name}-${idfile}.resquiggle.failed.tsv") into failed_resquiggles
    set folder_name, file ("${folder_name}-${idfile}"), file (".${folder_name}-${idfile}.RawGenomeCorrected_*.index") into fast5_folder_4_mods

    script:
	// conversion command if input is RNA
    def infolder = "${folder_name}-${idfile}"
    def demulti_cmd = "mkdir ${infolder}; cp ${fast5} ${infolder}"
    if (multi5 == "YES") {
        demulti_cmd = "echo 'Converting from multifast5 to single fast5'; mkdir ${infolder}; multi_to_single_fast5 -i ./ -s ${infolder}_tmp -t ${task.cpus}; rm ${infolder}_tmp/filename_mapping.txt; mv ${infolder}_tmp/*/*.fast5 ${infolder}"
    }    
    """ 
		#from multi to single
 	    ${demulti_cmd}
 	    # preprocessing adding 
        # resquiggling
        tombo resquiggle ${infolder}  ${reference} --rna --processes ${task.cpus} --overwrite --failed-reads-filename ${folder_name}-${idfile}.resquiggle.failed.tsv 
    """
}


fast5_groups_4_mods = fast5_folder_4_mods.groupTuple().into{
   fast5_groups_A_mods; fast5_groups_B_mods; fast5_groups_C_mods
}

/*
* Perform preprocessing for Epinano
*/
process index_reference {

	when:
	modfinder == "epinano"
    file(reference)

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

process epinano_get_variants {
    tag {"${modfinder}-${folder_name}"}  
    label 'big_mem_cpus'

	when:
	modfinder == "epinano"
	
    input:
    set folder_name, file(bam_file) from aln_4_epinano
    set file(reference), file(dict_index), file(faiidx) from indexes

    output:
    set folder_name, file("${folder_name}.tsv") into variants_for_splitting
   
    script:
	"""
	samtools view -h ${bam_file} -F 256 | \$SAM2TSV -r ${reference} >${folder_name}.tsv
	"""
}

process epinano_split_variants {
    tag {"${modfinder}-${folder_name}"}  

	when:
	modfinder == "epinano"
	
    input:
    set folder_name, file(variant) from variants_for_splitting

    output:
    set folder_name, file("${folder_name}_*.tsv") into variant_chunks
   
    script:
	"""
	split_tsv_file.py ${variant} ${params.granularity}
	"""	
}

/*
process epinano_get_per_read_variants {
    tag {"${modfinder}-${folder_name}"}  
    label 'big_mem_cpus'

	when:
	modfinder == "epinano"
	
    input:
    set folder_name, file(variant) from variant_chunks

    output:
    set folder_name, file("*.per_read.var.csv") into per_read_variants
   
    script:
	"""
	per_read_var.py ${variant} > ${}.per_read.var.csv
	"""
}

/*
process epinano_get_per_site_variants {
    tag {"${modfinder}-${folder_name}"}  
    label 'big_mem_cpus'

	when:
	modfinder == "epinano"
	
    input:
    set folder_name, file(variant) from variants_2_per_site

    output:
    set folder_name, file("${folder_name}.per_read.var.csv") into per_site_variants
    //set folder_name, file("${folder_name}.per_site.var.sliding.win.csv") into slide_per_site variants
   
    script:
	"""
	per_site_var.py ${folder_name}.tsv > ${folder_name}.per_site.var.csv
	#slide_per_site_var.py ${folder_name}.per_site.var.csv > ${folder_name}.per_site.var.sliding.win.csv
	"""
}

process epinano_get_events_from_fast5 {
    tag {"${modfinder}-${folder_name}"}  
    label 'big_mem_cpus'

	when:
	modfinder == "epinano2"
	
    input:
    set folder_name, file(fast5_batch) from fast5_for_epinano

    output:
    set folder_name, file("${fast5_batch}.event.tbl") into event_tbls
   
    script:
	"""
	fast5ToEventTbl.py ${fast5_batch} > ${fast5_batch}.event.tbl
	"""
}
	



/*
* detect modification 
*/

process getModifications {
    label 'big_mem_cpus'
	publishDir outputModifs, pattern: "*.significant_regions.fasta",  mode: 'copy'
        
    input:
    file(reference)
    set folder_name_A, file(fast5_folderA), file(indexA) from fast5_groups_A_mods.first()
    set folder_name_B, file(fast5_folderB), file(indexB) from fast5_groups_B_mods.last()
    
    output:
    file ("${folder_sample_name}_${folder_ctrl_name}.significant_regions.fasta") into sign_regions

    script:
    if (modfinder == "tombo") {
        def reference_cmd = unzipFile(reference, "reference.fa")
 	    """
 	    ${reference_cmd}
 	    mkdir ${folder_name_A} ${folder_name_B}
 	    mv ${fast5_folderA} ${folder_name_A}
 	    mv ${indexA} ${folder_name_A}
 	    mv ${fast5_folderB} ${folder_name_B}
 	    mv ${indexB} ${folder_name_B}
 	    tombo detect_modifications model_sample_compare --fast5-basedirs ${folder_sample_name}/* --control-fast5-basedirs ${folder_ctrl_name}/* --statistics-file-basename ${folder_sample_name}_${folder_ctrl_name}_model_sample_compare --rna --per-read-statistics-basename ${folder_sample_name}_${folder_ctrl_name}_per-read-statistics --processes ${task.cpus}
        tombo text_output signif_sequence_context --statistics-filename ${folder_sample_name}_${folder_ctrl_name}_model_sample_compare.tombo.stats  --genome-fasta reference.fa --fast5-basedirs ${folder_sample_name} --sequences-filename ${folder_sample_name}_${folder_ctrl_name}.significant_regions.fasta 
        rm reference.fa
        """
   } else if (modfinder == "epinano"){
		"""
		    
		"""	
	} else {
        """
 		echo "nothing to do!"
        """
   }
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
		R --slave -e "library(tailfindr); find_tails(fast5_dir = './' , save_dir = 'output', csv_filename = \'${folder_name}_findr.csv\', num_cores = ${task.cpus})"
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
		nanopolish polya -r ${fastq_file} -g ${reference_cmd} -t ${task.cpus} -b ${bam_file} > ${folder_name}.polya.estimation.tsv
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

// Get the channel for modfinder
def getModfinderChannel(batches, granularity, foldername, ismulti) {
	def number=0
	batches.collate(granularity).map { 
	    number++
		[foldername, number, ismulti, it]
	}.set{ fast5_4_modfinder }
	return fast5_4_modfinder
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

def unzipFile(zipfile, filename) {
    cmd = "ln -s ${zipfile} ${filename}"
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