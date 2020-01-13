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
BIOCORE@CRG Preprocessing of Nanopore direct RNA - N F  ~  version ${version}
====================================================

kit                       : ${params.kit}
flowcell                  : ${params.flowcell}
fast5                     : ${params.fast5}
reference                 : ${params.reference}
annotation                : ${params.annotation}

ref_type                  : ${params.ref_type}
seq_type                  : ${params.seq_type}

output                    : ${params.output}
qualityqc                 : ${params.qualityqc}
granularity               : ${params.granularity}

basecaller                : ${params.basecaller}
basecaller_opt            : ${params.basecaller_opt}
GPU                       : ${params.GPU}
demultiplexing            : ${params.demultiplexing} 
demultiplexing_opt        : ${params.demultiplexing_opt} 

filter                    : ${params.filter}
filter_opt                : ${params.filter_opt}
mapper                    : ${params.mapper}
mapper_opt                : ${params.mapper_opt}
map_type                  : ${params.map_type}

counter                   : ${ params.counter}
counter_opt               : ${ params.counter_opt}

email                     : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"
if (params.granularity == "") params.granularity = 1000000000

// check multi5 and GPU usage. GPU maybe can be removed as param if there is a way to detect it
if (params.GPU != "ON" && params.GPU != "OFF") exit 1, "Please specify ON or OFF in GPU processors are available"

// check sequence type parameter
//if (params.seq_type != "gDNA" && params.seq_type != "cDNA" && params.seq_type != "RNA") exit 1, "Please specify the sequence type as RNA, cDNA or gDNA"
//if (params.seq_type == "gDNA") { 
//    log.info "seqType is gDNA: map_type is set to 'unspliced'"
//	params.map_type = "unspliced"
//} else {
//    log.info "seqType is ${params.seq_type}: map_type is ${params.map_type}"
	if (params.map_type != "unspliced" && params.map_type != "spliced") exit 1, "Mapping type NOT supported! Please choose either 'spliced' or 'unspliced'"
//}

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
config_report = file("$baseDir/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("$baseDir/../docs/logo_small.png")
deeplexicon_folder = file("$baseDir/deeplexicon/")


basecaller   		= params.basecaller
basecaller_opt  	= params.basecaller_opt
demultiplexer 		= params.demultiplexing
demultiplexer_opt   = params.demultiplexing_opt
mapper      		= params.mapper
mapper_opt   		= params.mapper_opt
counter_opt   		= params.counter_opt 


// Output folders
outputFastq    = "${params.output}/fastq_files"
outputFast5    = "${params.output}/fast5_files"
outputQual     = "${params.output}/QC_files"
outputMultiQC  = "${params.output}/report"
outputMapping  = "${params.output}/alignment"
outputCounts   = "${params.output}/counts"
outputAssigned = "${params.output}/assigned"
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
    .into {fast5_4_name; fast5_4_testing; fast5_4_granularity}

/*
 * Get the name from the folder
 */
folder_info = params.fast5.tokenize("/")
folder_name = folder_info[-2]

// Check config file for consistency
//if (demultiplexer == "guppy" && params.barcodekit == "") 
//	exit 1, "Demultiplexing with guppy needs the definition of the barcodekit parameter. Exiting"
//if (basecaller != "guppy" && demultiplexer == "guppy") 
//	exit 1, "Demultiplexing with guppy can be performed ONLY when the basecaller is guppy too. Exiting"
//if (basecaller == "guppy" && demultiplexer == "guppy") 
//	log.info "Performing basecalling and demultiplexing at the same time with Guppy."

/*
* This is default value in case guppy will be used for RNA demultiplexing
*/
params.barcodekit = ""
if (demultiplexer == "") {
	demultiplexer = "OFF"
}

if (demultiplexer != "OFF" && demultiplexer != "deeplexicon")
exit 1, "Demultiplexing of RNA can be performed only with deeplexicon. Current value is ${demultiplexer}"

if (params.GPU == "YES" && basecaller != "guppy")
	exit 1, "GPU can be used only with GUPPY basecaller!"

if (params.ref_type == "genome") {
	annotation = file(params.annotation)
	if( !annotation.exists() ) exit 1, "Missing annotation file: ${params.annotation}! This is mandatory when ref_type = 'genome'"
}


process testInput {
    tag {"${fast5}"}  
            
    input:
    file(fast5) from fast5_4_testing.first()

    output:
    stdout into multi5_type

    script:
    """
    fast5_type.py ${fast5}
    """
}

multi5_type.map {  it.trim().toInteger() }.into{multi5_type_for_msg; multi5_type_for_bc; multi5_type_for_granularity; multi5_type_for_demultiplexing}
multi5_type_for_msg.map{it == 0 ? "Single Fast5 files detected!": "MultiFast5 files detected!" }.println()

// if you are using GPU analyse the whole dataset, otherwise make batch of 4,000 sequences if they are single fast5
// or single batches of multi fast5 sequences 
multi5_type_for_granularity.merge(fast5_4_granularity.collect()).map{
	(params.GPU == "YES" ? params.granularity  : (it[0] == 0 ? it[1..-1].collate(4000) : it[1..-1].collate(1)) )
}.flatMap().set{fast5_batches}

// create a map id batch -> list of files
def num_batch = -1
fast5_batches.map { 
    num_batch++
    [num_batch, it]
}.into{ fast5_4_basecall; fast5_4_demulti}

/*
*  Perform base calling using albacore or guppy on raw fas5 files
*/

process baseCalling {
    tag {"${basecaller}-${folder_name}-${idfile}"}  

    label (params.GPU == "ON" ? 'basecall_gpus': 'basecall_cpus')
	
	// move basecalled output fast5 files 
	if(demultiplexer != "deeplexicon") {
		publishDir outputFast5, pattern: "*_out/workspace/*.fast5",  mode: 'move', saveAs: { file -> "${file.split('\\/')[-1]}"  }
    }
    
    input:
    set idfile, file(fast5) from fast5_4_basecall
    val (multi5) from multi5_type_for_bc
    
    output:
    file ("${idfile}_out/workspace/*.fast5") optional true into fast5_files_for_demultiplexing
    set idfile, file ("${idfile}.*.gz") into fastq_files_for_demultiplexing, demulti_log
    file ("${idfile}_out/sequencing_summary.txt") into seq_summaries optional true 

    script:
	// conversion command if input is RNA - have to check if this is really needed
	def RNA_conv_cmd = ""
 	def demulti_cmd = ""
    def infolder = "./"
	if (params.seq_type == "RNA") {	RNA_conv_cmd = " | awk '{if (NR%4==2) gsub(\"U\",\"T\"); print}' " }   
    if (basecaller == "albacore") {
        // in case input files are multi fast5 convert them in single fast5 since albacore is not able to deal with multi fast5
        if (multi5 == 1) {
        	demulti_cmd = "mkdir demulti_tmp; mkdir demulti; multi_to_single_fast5 -i ${infolder} -s demulti_tmp -t ${task.cpus}; mv demulti_tmp/*/*.fast5 demulti; rm -fr demulti_tmp"
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
		mkdir single_basecallings
		mv ${idfile}_out/workspace/*/*.fast5 single_basecallings
		mkdir temp_multi
        single_to_multi_fast5 -i single_basecallings -s temp_multi -t ${task.cpus}
        mv temp_multi/batch_0.fast5 ./${idfile}_out/workspace/batch_${idfile}.fast5
		if [-d demulti]; then rm -fr demulti; fi
		rm -fr single_basecallings temp_multi
        """
   } else if (basecaller == "guppy"){
    	def multi_cmd = ""
        def gpu_cmd = ""
        def gpu_prefix = ""
        if (params.GPU == "ON") {
            gpu_prefix = 'export LD_LIBRARY_PATH="/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/.singularity.d/libs"'
         	gpu_cmd = '-x "cuda:0"'
        }
        // in case input files are single fast5 group them in multifast5 at the end
        if (multi5 == 0) {
            multi_cmd = "mkdir single_basecallings temp_multi; mv *_out/workspace/*.fast5 single_basecallings; single_to_multi_fast5 -i single_basecallings -s temp_multi -t ${task.cpus}; mv temp_multi/batch_0.fast5 ./${idfile}_out/workspace/batch_${idfile}.fast5; rm -fr temp_multi single_basecallings"
        }
        // Different command line in case guppy is also demultiplexer
        if (demultiplexer == "guppy") {
			"""
                ${gpu_prefix}
				guppy_basecaller ${gpu_cmd} ${basecaller_opt} ${demultiplexer_opt} --flowcell ${params.flowcell} --kit ${params.kit} --num_barcode_threads ${task.cpus} --barcode_kits ${params.barcodekit} --trim_barcodes  --fast5_out --input ${infolder} --save_path ./${idfile}_out --cpu_threads_per_caller 1  --num_callers ${task.cpus} 
				cd ${idfile}_out; 
				if [ -d barcode01 ]; then
					for d in barcode*; do echo \$d; cat \$d/*.fastq ${RNA_conv_cmd} > ../${idfile}.\$d.fastq; done;
				fi
				cat unclassified/*.fastq ${RNA_conv_cmd} > ../${idfile}.unclassified.fastq; cd ../
				for i in *.fastq; do gzip \$i; done
				${multi_cmd}
			"""
        }
        else  {
			"""
            ${gpu_prefix}
			guppy_basecaller ${gpu_cmd} --flowcell ${params.flowcell} --kit ${params.kit} --fast5_out ${basecaller_opt} --input ${infolder} --save_path ./${idfile}_out --cpu_threads_per_caller 1  --num_callers  ${task.cpus} 
			cat ${idfile}_out/*.fastq ${RNA_conv_cmd} >> ${idfile}.fastq
			rm ${idfile}_out/*.fastq
			gzip ${idfile}.fastq
			${multi_cmd}
			"""	
		}	
   } else if (demultiplexer == "OFF") {		
		"""
 		fast5_to_fastq.py ${fast5}; mv *.fastq ${idfile}.fastq; gzip ${idfile}.fastq
        """
   }
}

//demulti_log.println()

/*
*  Perform demultiplexing (optional) using porechop on basecalled reads
*/
if(demultiplexer == "deeplexicon") {
	process demultiplexing_with_deeplexicon {
		label 'basecall_cpus'
   	    tag {"${demultiplexer}-${idfile}"}  
				
		input:
    	set idfile, file(fast5) from fast5_4_demulti
    	val (multi5) from multi5_type_for_bc
        file(deeplexicon_folder)
        
		output:
		set idfile, file ("${idfile}_demux.tsv") into demux_for_fastq_extraction
		file ("${idfile}_demux.tsv") into demux_for_fast5_extraction

		script:
		def model = ''
		def deep_option = 'multi'
		if (multi5 == 0){
			deep_option = 'single'
		}
		"""
		    ln -s ${deeplexicon_folder}/* .
            deeplexicon.py -p ./ ${demultiplexer_opt} -f ${deep_option} -b 4000 -v > ${idfile}_demux.tsv
 		"""
	} 
	
	process extracting_demultiplexed_fastq {
		label 'basecall_cpus'
   	    tag {"${demultiplexer}"}  
				
		input:
    	set idfile, file(demux), file(fastq) from demux_for_fastq_extraction.join(fastq_files_for_demultiplexing)
        
		output:
		set idfile, file ("*.fastq.gz") into fastq_for_filtering

		script:
		"""
            extract_sequence_from_fastq.py ${demux} ${fastq}
			for i in *.fastq; do gzip \$i; done
 		"""
	} 
	
	process extracting_demultiplexed_fast5 {

		label 'basecall_cpus'
   	    tag {"${demultiplexer}-${idfile}"}  
		publishDir outputFast5,  mode: 'copy'
				
		input:
    	file("demux_*") from demux_for_fast5_extraction.collect()
    	file("*") from fast5_files_for_demultiplexing.collect()

		output:
        file("dem_*")
        
		script:
		"""
		cat demux_* | grep -v ReadID >> dem.files
		awk '{print \$2 > \$3".list" }' dem.files
		for i in *.list; do mkdir dem_`basename \$i .list`; done
		for i in *.list; do fast5_subset --input ./ --save_path dem_`basename \$i .list`/ --read_id_list \$i --batch_size 4000 -t ${task.cpus}; done 
 		rm dem_*/filename_mapping.txt
 		"""
	} 
	
	
} else {
	fastq_files_for_demultiplexing.set{ fastq_for_filtering}
}

/*
*  Perform filtering (optional) using nanofilt on fastq files
*/
if (params.filter == "nanofilt") {
	process filtering {
		label 'big_cpus'
   	    tag {"${params.filter}-${fastq_file}".replace('.fastq.gz', '')}  
				
		input:
		set idfile, file(fastq_file) from fastq_for_filtering.transpose()

		output:
		set idfile, file("*-filt.fastq.gz") into fastq_for_next_step

		script:
		output = "${fastq_file}".replace(".fastq.gz", "-filt.fastq.gz")
		"""
			zcat ${fastq_file} | NanoFilt ${params.filter_opt} | gzip > ${output}
		"""
	} 
} else {
	fastq_for_filtering.transpose().set{fastq_for_next_step}
}


fastq_for_next_step.map{
	filepath=it[1]
    if (demultiplexer != "OFF") {
        fileparts = filepath.getName().tokenize(".")
 		["${folder_name}.${fileparts[-3]}", filepath]
	} else {
		["${folder_name}", filepath]
	}
}.groupTuple().set{fastq_files_for_grouping}

/*
*  Concatenate FastQ files
*/
process concatenateFastQFiles {
    tag {idfile} 

	publishDir outputFastq, pattern: "*.fq.gz",  mode: 'copy'

    input:
    set idfile, file(fastq_files) from fastq_files_for_grouping

    output:
    set idfile, file("${idfile}.fq.gz") into fastq_files_for_fastqc, fastq_files_for_mapping

    script:
    """
		cat *.fastq.gz >> ${idfile}.fq.gz
    """
}

/*
*  Perform QC on fast5 files
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
*  Perform fastQC on fastq files
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
*  Perform counting (optional)
*/

if ( params.counter == "YES") {
	process counting {
		tag {"${idfile}"}  
		publishDir outputCounts, pattern: "*.count", mode: 'copy'
		publishDir outputAssigned, pattern: "*.assigned", mode: 'copy'

		input:
		set idfile, file(bamfile) from aligned_reads_for_counts

		output:
		file("${idfile}.count") into read_counts
		file("${idfile}.stats") optional true into count_stats
		file("${idfile}.assigned") optional true
		script:    
		if (params.ref_type == "transcriptome") {
			"""
			NanoCount -i ${bamfile} -o ${idfile}.count ${counter_opt};
	awk '{sum+=\$3}END{print FILENAME"\t"sum}' ${idfile}.count |sed s@.count@@g > ${idfile}.stats
	samtools view -F 256 ${bamfile} |cut -f 1,3 > ${idfile}.assigned
			"""
		} else if (params.ref_type == "genome") {
			def anno = unzipBash(annotation) 
			"""
			samtools view ${bamfile} |htseq-count -f sam - ${anno} -o ${idfile}.sam > ${idfile}.count
			awk '{gsub(/XF:Z:/,"",\$NF); print \$1"\t"\$NF}' ${idfile}.sam |grep -v '__' > ${idfile}.assigned
			rm ${idfile}.sam
			"""		
		}
	}

	/*
	*  Join alnQC 
	*/
	process joinCountQCs {
   
		input:
		file "*" from count_stats.collect()

		output:
		file("counts_mqc.txt") into count_repo_for_multiQC
	
		script:
		"""
	   echo '# id: NanoCount
	# plot_type: \'table\'
	# section_name: Read counts 
	File name	\'Counts\' ' > counts_mqc.txt 
		cat *.stats  >> counts_mqc.txt 
		"""
	}
} else {
   read_counts = Channel.empty()
   count_repo_for_multiQC = Channel.empty()
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
    publishDir outputQual, pattern: "*_plot/*", mode: 'copy'
    label 'big_cpus'
    errorStrategy 'ignore'
    tag {bamid}  
   
    input:
    set bamid, file(bamfile) from aligned_reads_for_QC2

    output:
    file("*_plot/*") optional true
    file("${bamid}_stats_mqc.png") optional true into qc2_for_multiqc
    
    script:
    """
    NanoPlot --bam ${bamfile} -o ${bamid}_plot --maxlength 5000 -t ${task.cpus}
    mkdir tmp_dir
    cp ${bamid}_plot/PercentIdentityvsAverageBaseQuality_kde.png tmp_dir
    cp ${bamid}_plot/LengthvsQualityScatterPlot_dot.png tmp_dir
    cp ${bamid}_plot/HistogramReadlength.png tmp_dir 
    cp ${bamid}_plot/Weighted_HistogramReadlength.png tmp_dir
    gm montage tmp_dir/*.png -tile 2x2 -geometry 800x800 ${bamid}_stats_mqc.png
    rm -fr tmp_dir
    """
}

QC_folders.mix(fastqc_for_multiqc,qc2_for_multiqc,read_counts,count_repo_for_multiQC,alnQC_for_multiQC).set{files_for_report}

/*
*  Perform multiQC report
*/
process multiQC {
    publishDir outputMultiQC, mode: 'copy'
   
    input:
    file(logo)
    file(config_report)
    file("*") from files_for_report.collect()
    
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

// make named pipe 
def unzipBash(filename) { 
    cmd = filename.toString()
    if (cmd[-3..-1] == ".gz") {
    	cmd = "<(zcat ${filename})"
    }
    return cmd
}

