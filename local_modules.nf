// LOCAL MODULES FOR MOP

params.LABEL = ""
params.OUTPUT = ""
params.saveSpace = "NO"

// MODULES 
// MOP_PREPROCESS
process extracting_demultiplexed_fastq {
    label (params.LABEL)
    tag "${ idfile }"
			
	input:
	tuple val(idfile), path(demux), path(fastq) 
	
	
	output:
	tuple val(idfile), path ("*.fastq.gz")

	script:
	"""
		extract_sequence_from_fastq.py ${demux} ${fastq}
		for i in *.fastq; do gzip \$i; done
	"""
}

process preparing_demultiplexing_fast5_deeplexicon {

    label (params.LABEL)
    tag "${ idfile }"
		
	input:
	tuple val(idfile), path("demux_*")

	output:
	tuple val(idfile), path("*.list")

	
	script:
	"""
	cat demux_* | grep -v ReadID >> dem.files
	awk '{print \$2 > \$3".list" }' dem.files
	"""
}

process extracting_demultiplexed_fast5_deeplexicon {
    label (params.LABEL)
	container 'lpryszcz/deeplexicon:1.2.0'
    tag "${ idfile } on ${ idlist }"
    if (params.saveSpace == "YES") publishDir(params.OUTPUTF5, mode:'move', pattern: '*-*') 
    else publishDir(params.OUTPUTF5, mode:'copy', pattern: '*-*')    

    publishDir(params.OUTPUTST, mode:'copy', pattern: 'summaries/*_final_summary.stats', saveAs: { file -> "${file.split('\\/')[-1]}" })    

		
	input:
	tuple val(idfile), path(idlist), file("*")

	output:
	path("${idfile}-*"), type: "dir", emit: dem_fast5
	path("summaries/*_final_summary.stats"), emit: dem_summaries
	
	script:
	"""
	mkdir ${idfile}---`basename ${idlist} .list`; fast5_subset --input ./ --save_path ${idfile}---`basename ${idlist} .list`/ --read_id_list ${idlist} --batch_size 4000 -c vbz -t ${task.cpus}
	mkdir summaries
	for i in */filename_mapping.txt; do awk 'BEGIN{print "filename\tread_id"}{print \$2"\t"\$1}' \$i > `echo \$i | awk -F"/" '{print "summaries/"\$1"_final_summary.stats"}'`; done
	rm */filename_mapping.txt;
	"""
} 

process extracting_demultiplexed_fast5_guppy {
    tag "${ idfile }"
    label (params.LABEL)
    if (params.saveSpace == "YES") publishDir(params.OUTPUT, mode:'move') 
    else publishDir(params.OUTPUT, mode:'copy')    

    container "quay.io/biocontainers/ont-fast5-api:4.0.0--pyhdfd78af_0"

             
	input:
	tuple val(idfile), path("summaries_*"), file("*")
    
	output:
	path("${idfile}-*")

    script:
    """
      if [ -f "summaries_" ]; then
	  ln -s summaries_ final_summary.stats
	  else 
		  head -n 1 summaries_1 > final_summary.stats
	      for i in summaries_*; do grep -v "filename" \$i | awk -F"\t" -v id=${idfile}  '{OFS="\t"; \$19 = id"---"\$19; print \$0}'  >> final_summary.stats; done
	  fi

		demux_fast5 -c vbz -t ${task.cpus} --input ./ --save_path ./ --summary_file final_summary.stats 
		rm -fr barcode_arrangement
    """
}


/*
*  Clean files
*/
process cleanFile {
    tag "${id}"
    label (params.LABEL)
    
    input:
    tuple val(id), path(file_to_remove)
    val(file_to_wait1)
	val(extension)

	when: params.saveSpace == "YES"
    
    script:
    """
		for i in *${extension}; do rm \$(readlink -f \$i); done
    """
}


/*
*  Concatenate FastQ files
*/
process concatenateFastQFiles {
    tag "${idfile}"
    label (params.LABEL)
    publishDir(params.OUTPUT, mode:'copy') 

    input:
    tuple val(idfile), path(demultifq)

    output:
    tuple val(idfile), path("${idfile}.fq.gz") 
    

    script:
    """
		cat *.fastq.gz  >> ${idfile}.fq.gz
    """
}

/*
*  Perform QC on fast5 files
*/

process MinIONQC {
    tag "${folder_name}"
    label (params.LABEL)
    container 'biocorecrg/mopprepr:0.7'
    errorStrategy 'ignore'
    if (params.OUTPUT != "") publishDir(params.OUTPUT, mode:'copy', pattern: '*.stats') 

    
    input:
    tuple val(folder_name), path("summaries_*") 

    output:
    tuple val(folder_name), path ("${folder_name}_QC"), emit: QC_folder 
    tuple val(folder_name), path ("*_summary.stats"), emit: stats 

    script:
    """
      if [ -f "summaries_" ]; then
	  ln -s summaries_ ${folder_name}_final_summary.stats
	  else 
		  head -n 1 summaries_1 > ${folder_name}_final_summary.stats
	      for i in summaries_*; do grep -v "filename" \$i >> ${folder_name}_final_summary.stats; done
	  fi
      MinIONQC.R -i ${folder_name}_final_summary.stats -o ${folder_name}_QC -q ${params.qualityqc} -p ${task.cpus}
    """
}

/*
*  Perform bam2stats QC 
*/
process bam2stats {
    label (params.LABEL)
    tag "${id}" 
   
    input:
    tuple val(id), path(bamfile)

    output:
    tuple val(id), path ("${id}.stat")
    
    script:
    """
    bam2stats.py ${bamfile} > ${id}.stat
    """
}

/*
*
*/

process AssignReads {
    tag "${id}"
    publishDir(params.OUTPUT, mode:'copy') 
    label (params.LABEL)

    input:
    tuple val(id), path(input)
	val(tool)

    output:
    tuple val(id), path ("${id}.assigned")
    
    script:
    if (tool == "nanocount")
	    """
		samtools view -F 256 ${input} |cut -f 1,3 > ${id}.assigned
    	"""
    else if(tool == "htseq")
    	"""
			samtools view ${input} | awk '{gsub(/XF:Z:/,"",\$NF); print \$1"\t"\$NF}' |grep -v '__' > ${id}.assigned
    	"""
    else 
        error "Invalid alignment mode: ${tool}"
}

/*
*
*/

process countStats {
    tag "${id}"
    label (params.LABEL)
   
    input:
    tuple val(id), path(input)

    output:
    tuple val(id), path ("${id}.count.stats")
    
    script:
	"""
		wc -l ${input} |sed s@.assigned@@g | awk '{print \$2"\t"\$1}' > ${id}.count.stats
    """
}

/*
*  Join AlnStats 
*/
process joinAlnStats {
    label (params.LABEL)
    tag "joining aln stats"
 
    input:
    file "alnqc_*" 

    output:
    path("alnQC_mqc.txt") 
    
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
*  Join Count Stats 
*/
process joinCountStats {
    tag "joining count stats"
    label (params.LABEL)
  
    input:
    file "stats_*" 

	output:
	path("counts_mqc.txt")
	
	script:
	"""
	echo '# id: Assigned reads
	# plot_type: \'table\'
	# section_name: Assigned counts 
	File name	\'Counts\' ' > counts_mqc.txt 
		cat stats_*  >> counts_mqc.txt 
		"""
} 

 process bam2Cram {
    tag "${idfile}"  
    
    publishDir(params.OUTPUT, mode:'copy') 
    label (params.LABEL)

    input:
    path(reference)
    val(subsampling_val)
    tuple val(idfile), path(aln), path(index)
    
    output:
    file("*.sorted.cram*") optional true 
    
    script:
    def downcmd = ""
    def input = aln
    def perc = 0
    def clean = ""
 	if (subsampling_val != "NO") {
		perc = subsampling_val/100
		downcmd = "samtools view -@ ${task.cpus} -bs ${perc} ${aln} > subsample.bam"
		input = "subsample.bam"
		clean = "rm subsample.bam"
	}
 	"""
 		${downcmd}
		samtools view  -@ ${task.cpus} -C ${input} -T ${reference} -o ${idfile}.sorted.cram
		samtools index -@ ${task.cpus} ${idfile}.sorted.cram
		${clean}
    """
}

// COMMON TO ALL
process checkRef {
    tag "Checking ${ reference }"
    label (params.LABEL)
 
    input:
    path(reference)
    
    output:
    path("reference.fa")
    
    script:
	"""
	if [ `echo ${reference} | grep ".gz"` ]; then 
   		zcat ${reference} > reference.fa
	else 
        ln -s ${reference} reference.fa
	fi
	"""	
}

// MOP_MOD and MOP_TAIL
process splitReference {
    label (params.LABEL)
    container 'biocorecrg/mopmod:0.6.2'
    tag "Splitting of ${ reference }"

    input:
    path(reference)

    output:
    path("pieces*.fa")

    script:
    """
    faSplit about ${reference} 100000 pieces
    """
}

process splitBams {
    label (params.LABEL)
    container 'biocorecrg/mopmod:0.6.2'
    tag "Splitting of ${ bams } on ${ref_piece}"

    input:
    tuple val(combid), path(bams), path(ref_piece)

    output:
    tuple val(combid), path("${combid}_s.bam"), path("${combid}_s.bam.bai")

    script:
    """
		samtools faidx ${ref_piece} 
		awk '{OFS="	"}{print \$1, "1", \$2}' ${ref_piece}.fai > ${ref_piece}.bed
		samtools view -@ ${task.cpus} ${bams} -L ${ref_piece}.bed -S | samtools view -Sb -t ${ref_piece}.fai -@ ${task.cpus} -o ${combid}.bam
   		samtools sort -@ ${task.cpus} -o ${combid}_s.bam ${combid}.bam
		samtools index ${combid}_s.bam
		rm ${combid}.bam
    """

}

process indexReference {
    label (params.LABEL)
    container 'biocorecrg/mopmod:0.6'
    tag "Indexing ${ reference }"

    input:
    path(reference)
    
    output:
    tuple val("${reference.simpleName}"), path(reference), path("*.dict"), path ("*.fai")
    
    script:
	"""
	\$PICARD CreateSequenceDictionary R=${reference} O=${reference}.dict
	samtools faidx ${reference}
	"""
}

process joinEpinanoRes {
    label (params.LABEL)
    container 'biocorecrg/mopmod:0.6.2'
    tag "joining on ${id}"
    publishDir(params.OUTPUT, mode:'copy') 

    input:
    tuple val(id), path(epinanores)
    
    output:
    tuple val(id), path("*.plus_strand.per.site.csv.gz"), emit: plusepi 
    tuple val(id), path("*.plus_strand.per.site.csv.gz"), emit: minusepi 

    
    script:
	"""
	if compgen -G "*.plus_strand.per.site.csv.gz" > /dev/null; then
		zcat *pieces*.plus_strand.per.site.csv.gz | awk '!(NR>1 && /#Ref/)' | gzip >>  ${id}.plus_strand.per.site.csv.gz
	fi
	if compgen -G "*.minus_strand.per.site.csv.gz" > /dev/null; then
		zcat *pieces*.minus_strand.per.site.csv.gz | awk '!(NR>1 && /#Ref/)' | gzip >>  ${id}.minus_strand.per.site.csv.gz
	fi	
	"""
}


/*
* 
*/

/*
* CALC MEAN PER POSITION
*/

process mean_per_pos {

    container 'biocorecrg/mopmod:0.7'
    label (params.LABEL)
    tag "${idsample}" 
	
    input:
    tuple val(idsample), path(event_align) 
    
    output:
    tuple val(idsample), path("*_perpos_median.parquet")


    script:
    
    """
	mean_per_pos.py -i ${event_align} -o `basename ${event_align} .fast5_event_align.tsv.gz`
	#gzip *_processed_perpos_median.tsv
    """
}

/*
* CONCAT MEAN PER POS
*/
process concat_mean_per_pos {

    container 'biocorecrg/mopmod:0.7'
    label (params.LABEL)
    tag "${idsample} on ${chr_file}" 
    	
    input:
    tuple val(idsample), path(event_align), path(chr_file)
    
    output:
    tuple val(idsample), path("${idsample}.gz")


    script:
    """
	Merging_processed_nanopolish_data.py -i *.parquet -o ${idsample}.gz -t ${task.cpus} -l ${chr_file}
    """
}


/*
* CONCAT CSV FILES 
*/
process concat_csv_files {

    container 'biocorecrg/mopmod:0.7'
    label (params.LABEL)
    tag "${idsample}" 
    
    publishDir(params.OUTPUT, mode:'copy') 
	
    input:
    tuple val(idsample), path("files_*")
    
    output:
    tuple val(idsample), path("${idsample}.csv.gz")


    script:
    """
	zcat files_* | awk '!(NR>1 && /reference_kmer/)'| pigz -p ${task.cpus} >> ${idsample}.csv.gz
    """
}

/*
*
*/

process callVariants {
    tag "${sampleID}" 
    container 'biocorecrg/mopmod:0.6'
    label (params.LABEL)
	
    input:
    tuple val(sampleID), path(alnfile), path(reference), path(dict_index), path(faiidx) 

    output:
    tuple val(sampleID), path("${sampleID}.tsv")
   
    script:
	"""
	samtools view -h ${alnfile} -F 256 | \$SAM2TSV -R ${reference} | cut -f 3 --complement  > ${sampleID}.tsv
	"""
}

process makeEpinanoPlots {
	publishDir params.OUTPUT, mode: 'copy'
	container "biocorecrg/mopnanotail:0.3"
    label (params.LABEL)

    tag {"${sampleIDA}--${sampleIDB} ${mode}"}  
	
    input:
    path(rscript)
    tuple val(sampleIDA), val(sampleIDB), path(per_site_varA), path(per_site_varB) 
    val(mode)
    
    output:
    path("*.pdf")
       
    script:
	"""
	Rscript --vanilla ${rscript} ${per_site_varA} ${sampleIDA} ${per_site_varB} ${sampleIDB} ${mode}  
	"""
}

process multiToSingleFast5 {
    container 'biocorecrg/mopmod:0.6'
    label (params.LABEL)

    tag "${idsample}"  
	
    input:
    tuple val(idsample), path(fast5)
    
    output:
    tuple val(idsample), path("${idsample}-single")
       
    script:
	"""
    mkdir ${idsample}-single;
    multi_to_single_fast5 -i ./ -s ./ -t ${task.cpus}; 
    rm ./filename_mapping.txt; 
    mv ./*/*.fast5 ${idsample}-single;
	"""
}

/*
*
*/
process bedGraphToWig {
    container 'biocorecrg/mopmod:0.6'
    tag "${idsample}"  
    errorStrategy 'ignore'
	
    input:
    path(chromsizes)
    tuple val(idsample), path(bedgraph)
    
    output:
    tuple val(idsample), path("*.bw")
       
    script:
    def ofname = "${bedgraph.baseName}.wig"
	"""
	bedgraph2wig.pl --bedgraph ${bedgraph} --wig ${ofname} --step 1 --compact
	wigToBigWig -clip ${ofname} ${chromsizes} ${bedgraph.baseName}.bw
	rm *.wig
	"""
}

/*
*/
process mergeTomboWigs {
    label (params.LABEL)
    tag "${combID}"  
	publishDir params.OUTPUT, pattern: "*_Tombo_Output.tsv.gz",  mode: 'copy'
	container "biocorecrg/mopmod:0.6"

   input:
    val(strand)
    tuple val(combID), path(coverage), path(covcontrol), path(statistic)

	output:
	path("*_Tombo_Output.tsv.gz") optional true 
	
	script:
	"""
	Merge_Tombo.py ${statistic} ${covcontrol} ${coverage} ${combID}.${strand}
	gzip *_Tombo_Output.tsv
	"""
}

/*
*/
process RNA2DNA {
    label (params.LABEL)
    tag "${id}"  

   input:
    tuple val(id), path(rnafqfile)

	output:
	tuple val(id), path("*_RNA.fq.gz") 
	
	script:
    def ofname = "${rnafqfile.baseName}_RNA.fq"

	"""
		RNA_to_DNA_fq.py -i ${rnafqfile} -o ${ofname}
		gzip ${ofname}
	"""
}


/*
*/
process wigToBigWig {
    label (params.LABEL)
    tag "${id}"  
	container "biocorecrg/mopmod:0.6"
    //errorStrategy 'ignore'

   input:
    path(chromsizes)
    tuple val(id), path(bedgraph)

	output:
	tuple val(id), path("*.bw") optional true 
	
	script:
    def ofname = "${bedgraph.baseName}.bw"

	"""
	size=`zcat ${bedgraph} | wc -l`;
   	if [[ \$size -gt 2 ]]
    then
		wigToBigWig -clip ${bedgraph} ${chromsizes} ${ofname}
    else
    	echo "empty wig"
    fi
	"""
}


// MOP_TAIL

process collect_tailfindr_results {
	publishDir params.OUTPUT, pattern: "*_findr.csv.gz",  mode: 'copy'
	tag "${ sampleID }"  
    label (params.LABEL)
	
	input:
	tuple val(sampleID), path("tailfin_*")
	
	output:
    tuple val(sampleID), path("${sampleID}.findr.len.gz"), emit: length 
    tuple val(sampleID), file ("*_findr.csv.gz"), emit: csv 

	script:
	"""
	zcat tailfin_* | awk '!(NR>1 && /tail_start/)' | gzip >>  ${sampleID}_findr.csv.gz
	zcat ${sampleID}_findr.csv | awk -F"," '{if (\$5!="" && \$1!="read_id") print \$1"\t"\$5}'| gzip > ${sampleID}.findr.len.gz
	"""

}

/*
*/
process join_nanotail_results {
	container "biocorecrg/mopnanotail:0.2"
    label (params.LABEL)
    tag "joining nanotail results"

    publishDir params.OUTPUT,  mode: 'copy'
	tag { sampleID }  
	
	input:
	tuple val(sampleID), path(nanopol), path(tailfindr), path(genes)
	file(joinScript)
	
	output:
	file("${sampleID}_*")
	
	script:
	"""
	Rscript --vanilla ${joinScript} ${tailfindr} ${nanopol} ${genes} ${sampleID}
	"""
	
}


/*
* Filter bams
*/
process filter_bam {
	tag "${ sampleID }"
    label (params.LABEL)
	
	input:
	file(reference)
	tuple val(sampleID), path(alignment)

	output:
	tuple val(sampleID), path("${sampleID}_filt.bam")

	script:
	"""
    #to keep only mapped reads and remove secondary alignments 
    samtools view -@ {task.cpus} -bF 260 ${alignment} > ${sampleID}_filt.bam 
	"""
}

// MOP CONSENSUS
process indexFasta {
    label (params.LABEL)

    tag "${reference}" 
	
    input:
    path(reference)
    
    output:
    stdout   
       
    script:
	"""
	samtools faidx ${reference}
	cut -f 1,2 ${reference}.fai
	"""
}

process getChromInfo {
    label (params.LABEL)

    tag "${reference}" 
	
    input:
    path(reference)
    
    output:
    path("chrom.sizes")   
       
    script:
	"""
	samtools faidx ${reference}
	cut -f 1,2 ${reference}.fai > chrom.sizes
	"""
}


process nanoConsensus {
	publishDir "${params.OUTPUT}/${sampleIDs}-${chrName}", mode: 'copy'
	container "biocorecrg/mop_consensus:0.1"
    label (params.LABEL)
    errorStrategy 'ignore'

    tag "${sampleIDs} on ${chrName}"  
	
    input:
    path(nanoConScript)
    path(nanoScripts)
    path(reference)
    tuple val(sampleIDs), path(Epi_Sample), path(Epi_IVT), path(NP_Sample), path(NP_IVT), path(Tombo), path(Nanocomp), val(chrName), val(chrStart), val(chrEnd)
    
    output:
    path("*")
       
    script:
	"""
	Rscript --vanilla ${nanoConScript} -Epi_Sample ${Epi_Sample} \
	 -Epi_IVT ${Epi_IVT} \
	 -NP_Sample ${NP_Sample} \
	 -NP_IVT ${NP_IVT} \
	 -Tombo ${Tombo} \
	 -Nanocomp ${Nanocomp} \
	 -chr ${chrName} \
	 -ini_pos ${chrStart} -fin_pos ${chrEnd} \
	 -output ${sampleIDs} \
	 -fasta ${reference} \
	 --nanocomp_stat GMM_logit_pvalue 
	"""
}


/*
* COMMON FUNCTIONS
*/

// Create a channel for tool options
def getParameters(pars_tools_file) {
	def pars_tools = file(pars_tools_file)
	if( !pars_tools.exists() ) exit 1, "Missing tools options config: '$pars_tools'"

	def progPars = [:]
	def allLines  = pars_tools.readLines()

	for( line : allLines ) {
    	def list = line.split("\t")
    	if (list.length <3) {
			 error "ERROR!!! Tool option file has to be tab separated\n" 
		}
    	if (!(list[0] =~ /#/ )) {
			progPars["${list[0]}--${list[1]}"] = list[2].replace("\"", "").replace('$baseDir', "${baseDir}").replace('${baseDir}', "${baseDir}")
    	}  
	}	
	return(progPars)
}

// Create a channel for tool options
def parseFinalSummary(final_summary_file) {
	def outstring = ""
	if (final_summary_file != "" ) {
		final_summary = file(final_summary_file)
		if( final_summary.exists() ) {
			def allLines  = final_summary.readLines()

			for( line : allLines ) {
				def list = line.split("=")
				if (list[0] == "protocol") {
					def vals = list[1].split(":")
					outstring = "--flowcell ${vals[1]} --kit ${vals[2]}"
				}  
			}	
		} else {
			log.info '***No configuration file found!!. You must specify kit and flowcell in the parameters!!***\n'
			} 
		} else {
			log.info '***No configuration file given!!. You must specify kit and flowcell in the parameters!!***\n'
		}
	return(outstring)
}


def reshapeDemuxSamples(inputChannel) {
	def reshapedChannel = inputChannel.map {
 		def ids = it[0].split("---")
 		def dems = it[0].split("\\.")
 			["${ids[0]}---${dems[-1]}", it[1]]
	 }
	 return(reshapedChannel)
}

def reshapeSamples(inputChannel) {
    def reshapedChannel = inputChannel.map{
		def ids = it[0].split("---")
		["${ids[0]}", it[1]]
	}
	return(reshapedChannel)
}

def mapIDPairs (ids, values) {
	def combs = ids.combine(values, by:0).map{
		[it[1], it[0], it[2]]
	}.combine(values, by:0).map{
		[it[1], it[0], it[2],  it[3]]
	}
	return(combs)
}

def checkTools(tool_names, tool_lists) {
	println "----------------------CHECK TOOLS -----------------------------"
	tool_names.each{ key, value -> 
		if (value == "NO" ) {
			println "> ${key} will be skipped"
		} else {
			def combid = "${key}--${value}".toString()
			if (tool_lists.containsKey(combid)) {
				println "${key} : ${value}"	
			} else {
				println "ERROR ################################################################"
				println "${value} is not a valid program for ${key}"
				println "ERROR ################################################################"
				println "Exiting ..."
				System.exit(0)
			}
		}
	}
	println "--------------------------------------------------------------"
}
