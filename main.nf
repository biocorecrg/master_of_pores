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
BIOCORE@CRG Preprocessing of Nanopore data (gDNA, cDNA or RNA) - N F  ~  version ${version}
====================================================

kit                       : ${params.kit}
flowcell                  : ${params.flowcell}
fast5                     : ${params.fast5}
multi5                    : ${params.multi5}
reference                 : ${params.reference}

seqtype                   : ${params.seqtype}
output                    : ${params.output}
granularity               : ${params.granularity}
qualityqc                 : ${params.qualityqc}

basecaller                : ${params.basecaller}
basecaller_opt            : ${params.basecaller_opt}
GPU                       : ${params.GPU}
demultiplexing            : ${params.demultiplexing} 
demultiplexing_opt        : ${params.demultiplexing_opt} 
barcodekit                : ${params.barcodekit}
filter                    : ${params.filter}
filter_opt                : ${params.filter_opt}
mapper                    : ${params.mapper}
mapper_opt                : ${params.mapper_opt}
map_type                  : ${params.map_type}

email                     : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// check multi5 and GPU usage. GPU maybe can be removed as param if there is a way to detect it
if (params.multi5 != "YES" && params.multi5 != "NO") exit 1, "Please specify YES or NOT in case the data is in multi fast5 format or sigle fast5"
if (params.GPU != "ON" && params.GPU != "OFF") exit 1, "Please specify ON or OFF in GPU processors are available"

// check sequence type parameter
if (params.seqtype != "gDNA" && params.seqtype != "cDNA" && params.seqtype != "RNA") exit 1, "Please specify the sequence type as RNA, cDNA or gDNA"


if (params.seqtype == "gDNA") { 
    log.info "seqType is gDNA: map_type is set to 'unspliced'"
	params.map_type = "unspliced"
} else {
    log.info "seqType is ${params.seqtype}: map_type is ${params.map_type}"
	if (params.map_type != "unspliced" && params.map_type != "spliced") exit 1, "Mapping type NOT supported! Please choose either 'spliced' or 'unspliced'"
}

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
config_report = file("$baseDir/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("$baseDir/docs/logo_small.png")

basecaller   		= params.basecaller
basecaller_opt  	= params.basecaller_opt
demultiplexer 		= params.demultiplexing
demultiplexer_opt   = params.demultiplexing_opt
mapper      		= params.mapper
mapper_opt   		= params.mapper_opt

// Output folders
outputFastq    = "${params.output}/fastq_files"
outputFast5    = "${params.output}/fast5_files"
outputTar      = "${params.output}/fast5_tar_files"
outputQual     = "${params.output}/QC"
outputMultiQC  = "${params.output}/report"
outputMapping  = "${params.output}/aln"
outputCounts   = "${params.output}/counts"
outputReport   = file("${outputMultiQC}/multiqc_report.html")

/*
* move old multiQCreport
*/
if( outputReport.exists() ) {
  log.info "Moving old report to multiqc_report.html multiqc_report.html.old"
  outputReport.moveTo("${outputMultiQC}/multiqc_report.html.old")
}

/*
 * Creates the channels that emits fast5 files
 */
Channel
    .fromPath( params.fast5)                                             
    .ifEmpty { error "Cannot find any file matching: ${params.fast5}" }
    .into { fast5_4_to_batches; fast5_4_tar; fast5_4_name}

/*
 * Get the name from the folder
 */
folder_info = params.fast5.tokenize("/")
folder_name = folder_info[-2]


// Check config file for consistency
if (basecaller != "guppy" && demultiplexer == "guppy") 
	exit 1, "Demultiplexing with guppy can be performed ONLY when the basecaller is guppy too!!!\nERROR ~ Please change the config file!!"
if (basecaller == "guppy" && demultiplexer == "guppy") 
	log.info "Performing basecalling and demultiplexing at the same time with Guppy"	

def num_batch = 0
fast5_4_to_batches.collate(params.granularity).map { 
	num_batch++
	[num_batch, it]
}.set{ fast5_4_basecall }

/*
*  Perform base calling using albacore or guppy on raw fas5 files
*/

process baseCalling {
    tag {"${basecaller}-${folder_name}-${idfile}"}  

    label (params.GPU == "ON" ? 'basecall_gpus': 'basecall_cpus')
	
	publishDir outputFast5, pattern: "*_out/workspace/*.fast5",  mode: 'copy', saveAs: { file -> if (params.multi5 == "YES")  "${folder_name}/${file.split('\\/')[-1]}"  }
	publishDir outputTar, pattern: "*.tar",  mode: 'copy' 
            
    input:
    set idfile, file(fast5) from fast5_4_basecall

    output:
    file ("${idfile}_out/workspace/*.fast5") optional true
    file ("*.tar") optional true into fast5_files_for_tar
    file ("${idfile}.*.gz") into fastq_files_for_demultiplexing
    file ("${idfile}_out/sequencing_summary.txt") into seq_summaries

    script:
	// conversion command if input is RNA - have to check if this is really needed
	def RNA_conv_cmd = ""
	if (params.seqtype == "RNA") {	RNA_conv_cmd = " | awk '{if (NR%4==2) gsub(\"U\",\"T\"); print}' " }   
    def tar_cmd = "sleep 30; mkdir ${folder_name}_${idfile}_f5; mv *_out/workspace/*/*.fast5 ${folder_name}_${idfile}_f5; tar -cf ${folder_name}_${idfile}_f5.tar ${folder_name}_${idfile}_f5"
 	
    def infolder = "./"
    if (basecaller == "albacore") {
        def demulti_cmd = ""
        if (params.multi5 == "YES") {
        	demulti_cmd = "mkdir demulti; multi_to_single_fast5 -i ${infolder} -s demulti_tmp -t ${task.cpus}; rm demulti_tmp/filename_mapping.txt; mv demulti_tmp/*/*.fast5 demulti"
            infolder = "demulti"
        }
 	    """
 	    ${demulti_cmd}
 	    export PYTHONPATH=$baseDir/bin/albacore:\$PYTHONPATH
		read_fast5_basecaller.py ${basecaller_opt} --flowcell \"${params.flowcell}\" --kit \"${params.kit}\" --output_format fastq,fast5 \
			--worker_threads ${task.cpus} -s ./${idfile}_out --disable_filtering --input ${infolder};
		cat ${idfile}_out/workspace/*.fastq ${RNA_conv_cmd} >> ${idfile}.fastq
		rm ${idfile}_out/workspace/*.fastq
		gzip ${idfile}.fastq
        ${tar_cmd}
        """
   } else if (basecaller == "guppy"){
 		if (params.multi5 == "YES") { 
 			tar_cmd = ""
 		} 	
        def gpu_cmd = ""
        if (params.GPU == "ON") {
            
         	gpu_cmd = '-x "cuda:0"'
        }
        // Different command line in case guppy is also demultiplexer
        if (demultiplexer == "guppy") {
			"""
				guppy_basecaller ${gpu_cmd} ${basecaller_opt} ${demultiplexer_opt} --flowcell ${params.flowcell} --kit ${params.kit} --num_barcode_threads ${task.cpus} --barcode_kits ${params.barcodekit} --trim_barcodes  --fast5_out --input ${infolder} --save_path ./${idfile}_out --cpu_threads_per_caller 1  --num_callers ${task.cpus} 
				cd ${idfile}_out; 
				if [ -d barcode01 ]; then
					for d in barcode*; do echo \$d; cat \$d/*.fastq ${RNA_conv_cmd} > ../${idfile}.\$d.fastq; done;
				fi
				cat unclassified/*.fastq ${RNA_conv_cmd} > ../${idfile}.unclassified.fastq; cd ../
				for i in *.fastq; do gzip \$i; done
				${tar_cmd}
			"""
        }
        else {
			"""
			guppy_basecaller ${gpu_cmd} --flowcell ${params.flowcell} --kit ${params.kit} --fast5_out ${basecaller_opt} --input ${infolder} --save_path ./${idfile}_out --cpu_threads_per_caller 1  --num_callers  ${task.cpus} 
			cat ${idfile}_out/*.fastq ${RNA_conv_cmd} >> ${idfile}.fastq
			rm ${idfile}_out/*.fastq
			gzip ${idfile}.fastq
			${tar_cmd}
			"""	
		}	
   } else {
        """
 		echo "nothing to do!"
        """
   }
}


/*
*  Perform demultiplexing (optional) using porechop on basecalled reads
*/
if (demultiplexer == "porechop") {
	process demultiplexing_with_porechop {
		label 'basecall_cpus'
   	    tag {"${demultiplexer}-${fastq_file}"}  
				
		input:
		file(fastq_file) from fastq_files_for_demultiplexing

		output:
		file "*.fastq" into fastq_for_filtering

		script:
		"""
			porechop -i ${fastq_file} --threads ${task.cpus} -b ./ ${demultiplexer_opt}
		"""
		} 
	}
else{
	fastq_files_for_demultiplexing.set{ fastq_for_filtering}
}

/*
*  Perform filtering (optional) using nanofilt on fastq files
*/
if (params.filter == "nanofilt") {
	process filtering {
		label 'big_cpus'
   	    tag {"${params.filter}-${fastq_file}"}  
				
		input:
		file(fastq_file) from fastq_for_filtering

		output:
		file "*.filt.fastq.gz" into fastq_for_next_step

		script:
		output = "${fastq_file.getSimpleName()}.filt.fastq"
		"""
			zcat ${fastq_file} | NanoFilt ${params.filter_opt} | gzip > ${output}.gz
		"""
	} 
} else {
	fastq_for_filtering.set{fastq_for_next_step}
}



fastq_for_next_step.map{
	filepath=it
    fileparts=filepath.getSimpleName()
    if (demultiplexer != "") {
 		["${folder_name}.${fileparts}", filepath]
	} else {
		["${folder_name}", filepath]
	}
}.set{fastq_files_for_grouping}


/*
*  Concatenate FastQ files
*/
process concatenateFastQFiles {
    tag {idfile}  
	publishDir outputFastq, pattern: "*.fq.gz",  mode: 'copy'
    
    input:
    set idfile, file(fastq_files) from fastq_files_for_grouping.groupTuple()

    output:
    set idfile, file("${idfile}.fq.gz") into fastq_files_for_fastqc, fastq_files_for_mapping

    script:
    """
		cat *.fastq.gz >> ${idfile}.fq.gz
    """
}

//

/*
*  Perform tar on fast5 files
*/

process collectMultiFast5 {
    tag { folder_name }  
    publishDir outputTar, mode: 'copy'
    
    input:
    file(fast5) from fast5_files_for_tar.collect()

    when: 
    params.multi5 != "YES"

    output:
    file "${folder_name}.tar"

    script:
    """
       tar -hcvf ${folder_name}.tar ${fast5}
    """
}

/*
*  Perform QC on tar on fas5 files
*/

process QC {
    tag {folder_name}  
    label 'big_cpus'
    publishDir outputQual, mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    file("summaries_*") from seq_summaries.collect()

    output:
    file ("${folder_name}_QC") into QC_folders

    script:
    """
      if [ -f "summaries_" ]; then
	  ln -s summaries_ final_summary.stats
	  else 
		  head -n 1 summaries_1 > final_summary.stats
	      for i in summaries_*; do grep -v "filename" \$i >> final_summary.stats; done
	  fi
      MinIONQC.R -i final_summary.stats -o ${folder_name}_QC -q ${params.qualityqc} -p ${task.cpus}
    """
}

/*
*  Perform fastQC on tar on fastq files
*/

process fastQC {
    tag {idfile}  
    label 'big_cpus'

    publishDir outputQual, pattern: "*_fastqc.html", mode: 'copy'
   
    input:
    set idfile, file(fastq_file) from fastq_files_for_fastqc

    output:
    file ("*_fastqc.*") into fastqc_for_multiqc

    script:
    """
     fastqc ${fastq_file} -t ${task.cpus}
    """
}

/*
*  Perform mapping and sorting
*/
process mapping {
    tag {"${mapper}-${idfile}"}  
    publishDir outputMapping, mode: 'copy'
    label 'big_mem_cpus'

    input:
    file(reference)
    set idfile, file (fastq_file) from fastq_files_for_mapping
    
    output:
    set idfile, file("${idfile}.${mapper}.sorted.bam") optional true into aligned_reads, aligned_reads_for_QC, aligned_reads_for_QC2, aligned_reads_for_counts

    script:    
    if (mapper == "minimap2") {
	    def mappars = (params.map_type == "spliced") ? "-ax splice -k14" : "-ax map-ont"
	    mappars += " ${mapper_opt} "
 	    """
        minimap2 -t ${task.cpus} ${mappars} -uf ${reference} ${fastq_file} | samtools view -@ ${task.cpus} -F4 -hSb - > reads.mapped.bam
        samtools sort -@ ${task.cpus} -o ${idfile}.${mapper}.sorted.bam reads.mapped.bam
        rm reads.mapped.bam
        """
   }
   else if (mapper == "graphmap2"){
	    def mappars = (params.map_type == "spliced") ? "-x rnaseq" : ""
 	    mappars += " ${mapper_opt} "
        """
        graphmap2 align -t ${task.cpus} -r ${reference} ${mappars} -d ${fastq_file}  | samtools view -@ ${task.cpus} -F4 -hSb - > reads.mapped.bam
        samtools sort -@ ${task.cpus} -o ${idfile}.${mapper}.sorted.bam reads.mapped.bam
        rm reads.mapped.bam
        """
   } else {
        """
 		echo "nothing to do!"
        """
   }     
}

/*
*  Perform alnQC 
*/
process alnQC {
    tag {bamid}  
   
    input:
    set bamid, file(bamfile) from aligned_reads_for_QC

    output:
    file "${bamid}.stat" into single_alnQC_outs
    
    script:
    """
    bam2stats.py ${bamfile} > ${bamid}.stat
    """
}

/*
*  Join alnQC 
*/
process joinAlnQCs {
   
    input:
    file "alnqc_*" from single_alnQC_outs.collect()

    output:
    file("alnQC_mqc.txt") into alnQC_for_multiQC
    
    script:
    """
    echo '# id: alnQC
# plot_type: \'table\'
# section_name: \'Alignment QC\' ' > alnQC_mqc.txt
    cat alnqc_* | head -n 1| sed s@#@@g >> alnQC_mqc.txt
    cat alnqc_* | grep -v "#" >> alnQC_mqc.txt
    """
}

/*
*  Perform alnQC2 
*/

process alnQC2 {
    publishDir outputQual, pattern: "${bamid}_plot/", mode: 'copy'
    errorStrategy 'ignore'
    tag {bamid}  
   
    input:
    set bamid, file(bamfile) from aligned_reads_for_QC2

    output:
    file ("${bamid}_plot") optional true
    set file ("${bamid}_plot/PercentIdentityvsAverageBaseQuality_kde.png"), file ("${bamid}_plot/LengthvsQualityScatterPlot_dot.png"), file ("${bamid}_plot/HistogramReadlength.png"), file ("${bamid}_plot/Weighted_HistogramReadlength.png")  optional true into qc2_for_multiqc
    
    script:
    """
    NanoPlot --bam ${bamfile} -o ${bamid}_plot --title ${bamid} --maxlength 5000
    """
}


/*
*  Perform multiQC report
*/

process multiQC {
    publishDir outputMultiQC, mode: 'copy'
    
    input:
    file("*") from QC_folders.collect()
    file("*") from fastqc_for_multiqc.collect()
    set file("Percent_identity_vs_Average_base_quality_mqc.png"), file("Length_vs_Quality_dot_mqc.png"), file("Read_length_distribution_mqc.png"), file("Weighted_read_length_distribution_mqc.png") from qc2_for_multiqc
    file(alnQC_for_multiQC)
    file(config_report)
    file(logo)
    
    output:
    file("multiqc_report.html") into multiQC 
    
    script:
    """
     multiqc -c ${config_report} .
    """
}


if (params.email == "yourmail@yourdomain" || params.email == "") { 
    log.info 'Skipping the email\n'
}
else {
    log.info "Sending the email to ${params.email}\n"

    workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()

        sendMail(to: params.email, subject: "Master of Pore execution", body: msg,  attach: "${outputMultiQC}/multiqc_report.html")
    }
}

workflow.onComplete {
    println "Pipeline BIOCORE@CRG Master of Pore completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}


