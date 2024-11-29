#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '3.0'

params.help            = false
params.resume          = false

log.info """

╔╦╗╔═╗╔═╗  ╔═╗┬─┐┌─┐┌─┐┬─┐┌─┐┌─┐┌─┐┌─┐┌─┐
║║║║ ║╠═╝  ╠═╝├┬┘├┤ ├─┘├┬┘│ ││  ├┤ └─┐└─┐
╩ ╩╚═╝╩    ╩  ┴└─└─┘┴  ┴└─└─┘└─┘└─┘└─┘└─┘

====================================================
BIOCORE@CRG Master of Pores 3. Preprocessing - N F  ~  version ${version}
====================================================

conffile                  : ${params.conffile}

fast5                     : ${params.fast5}
fastq                     : ${params.fastq}

reference                 : ${params.reference}
annotation                : ${params.annotation}

granularity               : ${params.granularity}
ref_type                  : ${params.ref_type}

pars_tools                : ${params.pars_tools}
barcodes                  : ${params.barcodes}

output                    : ${params.output}

GPU                       : ${params.GPU}

basecalling               : ${params.basecalling}
demultiplexing            : ${params.demultiplexing}
demulti_fast5             : ${params.demulti_fast5}

filtering                 : ${params.filtering}
mapping                   : ${params.mapping}

counting                  : ${params.counting}
discovery                 : ${params.discovery}

cram_conv                 : ${params.cram_conv}
subsampling_cram          : ${params.subsampling_cram}

email                     : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// include functions, outdirs from other files
evaluate(new File("../outdirs.nf"))
//def local_modulesDir = "${projectDir}/../local"
def local_modules = file("${projectDir}/../local_modules.nf")
def subworkflowsDir = "${projectDir}/../BioNextflow/subworkflows"
def workflowsDir = "${projectDir}/../BioNextflow/workflows"
joinScript = file("${projectDir}/bin/join.r")

// get and check input files
if (params.mapping != "NO") {
    reference = file(params.reference)
    if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
} else {
    reference = ""
}

// INIZIALIZE MULTIQC REPORT
config_report = file("${projectDir}/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("${projectDir}/../img/logo_small.png")

Channel.from( config_report, logo ).set{multiqc_data}

outputReport   = file("${outputMultiQC}/multiqc_report.html")

if( outputReport.exists() ) {
  log.info "Moving old report to multiqc_report.html multiqc_report.html.old"
  outputReport.moveTo("${outputMultiQC}/multiqc_report.html.old")
}

// Get models

demux_models = ""

switch(params.demultiplexing) {
    case "deeplexicon":
        demux_models = "${projectDir}/deeplexicon_models/"
    break;
    case "seqtagger":
        demux_models = "${projectDir}/seqtagger_models/"
    break;
}

dorado_models = "${projectDir}/dorado_models/"

// check GPU usage.
if (params.GPU != "cuda11" && params.GPU != "cuda10" && params.GPU != "OFF" && params.GPU != "ON") exit 1, "Please specify cuda11, cuda10, ON or OFF if GPU processors are available. ON is legacy for cuda10"
def gpu = (params.GPU != 'OFF' ? 'ON' : 'OFF')
def cuda_cont = (params.GPU == 'cuda11' ? 'biocorecrg/mopbasecallc11:0.3' : 'biocorecrg/mopbasecall:0.3')


// CHECK INCOMPATIBILITIES AMONG PARAMETERS

if (params.ref_type == "genome") {
    if (params.annotation != "") {
        annotation = file(params.annotation)
        if( !annotation.exists() ) exit 1, "Missing annotation file: ${params.annotation}!"
    }
}

outmode = "copy"

include { final_message; notify_slack } from "${subworkflowsDir}/global_functions.nf"
include { checkInput; filterPerBarcodes; get_barcode_list; RNA2DNA; parseFinalSummary; checkTools; reshapeSamples; reshapeDemuxSamples; checkRef; getParameters; homogenizeVals } from "${local_modules}"

def demulti_fast5_opt = homogenizeVals(params.demulti_fast5)
def basecall_label = (params.GPU != 'OFF' ? 'basecall_gpus' : 'big_cpus')
def deeplexi_basecall_label = (params.GPU != 'OFF' ? 'demulti_gpus' : '')


def output_bc = (demulti_fast5_opt == 'ON' ? '' : outputFast5)
//def outputMinionQC = (demulti_fast5_opt == 'ON' ? '': outputQual)


def guppypars = ""
// GET PROGRAM PARS AND VERIFY
def tools = [:]
tools["basecalling"] = homogenizeVals(params.basecalling)
tools["demultiplexing"] = homogenizeVals(params.demultiplexing)
tools["mapping"] = homogenizeVals(params.mapping)
tools["filtering"] = homogenizeVals(params.filtering)
tools["counting"] = homogenizeVals(params.counting)
tools["discovery"] = homogenizeVals(params.discovery)

// Remove basecalling and demultiplexing in case of fastq input
if(params.fast5 == "" && params.fastq != "") {
    tools["basecalling"] = "NO"
    tools["demultiplexing"] = "NO"
} else {
    guppypars = parseFinalSummary(params.conffile)
    // Create a channel for tool options
    if (workflow.profile == "awsbatch") guppypars = guppypars + " --data_path /nextflow-bin/ont-guppy/data"
}

progPars = getParameters(params.pars_tools)
checkTools(tools, progPars)

// Create a channel for excluded ids
barcodes_to_include = get_barcode_list(params.barcodes)

def guppy_basecall_pars = guppypars + " " + progPars["basecalling--guppy"]

def basecaller_pars = ["guppy" : guppy_basecall_pars, "dorado" : progPars["basecalling--dorado"] ]
def demux_pars = ["guppy" : progPars["demultiplexing--guppy"] + " " + guppy_basecall_pars, "seqtagger":  progPars["demultiplexing--seqtagger"], "deeplexicon": progPars["demultiplexing--deeplexicon"] ]


// INCLUDE WORKFLOWS
include { BASECALL } from "${workflowsDir}/basecaller" addParams(gpu: gpu, output: output_bc, label: basecall_label, type:params.basecalling ,  extrapars: basecaller_pars[params.basecalling], models: dorado_models )
include { DEMULTIPLEX } from "${workflowsDir}/demultiplexer.nf" addParams(gpu: gpu, output: output_bc, label: basecall_label, type:params.demultiplexing , extrapars: demux_pars[params.demultiplexing], models: demux_models )
include { BASECALL_DEMULTIPLEX } from "${workflowsDir}/basecaller_demultiplexer.nf" addParams(gpu: gpu, output: output_bc, label: basecall_label, type:params.demultiplexing , extrapars: demux_pars[params.demultiplexing] )
include { DEMULTI_FAST5; DEMULTI_FAST5_FILTER } from "${subworkflowsDir}/misc/demulti_fast5" addParams(OUTPUT: outputFast5, OUTPUTST: outputQual, LABEL: 'big_cpus', TYPE: params.demultiplexing)


// INCLUDE MODULES
include { GET_VERSION as NANOFILT_VER; FILTER as NANOFILT_FILTER} from "${subworkflowsDir}/trimming/nanofilt" addParams(EXTRAPARS: progPars["filtering--nanofilt"])
include { GET_VERSION as NANOQ_VER; FILTER as NANOQ_FILTER} from "${subworkflowsDir}/trimming/nanoq" addParams(EXTRAPARS: progPars["filtering--nanoq"])
include { MAP as GRAPHMAP} from "${subworkflowsDir}/alignment/graphmap" addParams(EXTRAPARS: progPars["mapping--graphmap"], LABEL:'big_mem_cpus')
include { MAP as GRAPHMAP2} from "${subworkflowsDir}/alignment/graphmap2" addParams(EXTRAPARS: progPars["mapping--graphmap2"], LABEL:'big_mem_cpus')
include { MAP as MINIMAP2} from "${subworkflowsDir}/alignment/minimap2" addParams(EXTRAPARS: progPars["mapping--minimap2"], LABEL:'big_mem_cpus')
include { ALL as BWA} from "${subworkflowsDir}/alignment/bwa" addParams(EXTRAPARS: progPars["mapping--bwa"], LABEL:'big_mem_cpus')
include { GET_VERSION as BWA_VER} from "${subworkflowsDir}/alignment/bwa"
include { GET_VERSION as GRAPHMAP_VER} from "${subworkflowsDir}/alignment/graphmap"
include { GET_VERSION as GRAPHMAP2_VER} from "${subworkflowsDir}/alignment/graphmap2"
include { GET_VERSION as MINIMAP2_VER} from "${subworkflowsDir}/alignment/minimap2"
include { FASTQCP as FASTQC} from "${subworkflowsDir}/qc/fastqc" addParams(LABEL: 'big_cpus')
include { GET_VERSION as FASTQC_VER} from "${subworkflowsDir}/qc/fastqc"
include { SORT as SAMTOOLS_SORT } from "${subworkflowsDir}/misc/samtools" addParams(LABEL: 'big_cpus', OUTPUT:outputMapping)
include { INDEX as SAMTOOLS_INDEX } from "${subworkflowsDir}/misc/samtools" addParams(OUTPUT:outputMapping)
include { GET_VERSION as SAMTOOLS_VERSION; CAT as SAMTOOLS_CAT } from "${subworkflowsDir}/misc/samtools"
include { MOP_QC as NANOPLOT_QC } from "${subworkflowsDir}/qc/nanoplot" addParams(LABEL: 'big_cpus_retry')
include { GET_VERSION as NANOPLOT_VER } from "${subworkflowsDir}/qc/nanoplot"
include { GET_VERSION as NANOCOUNT_VER } from "${subworkflowsDir}/read_count/nanocount"
include { COUNT as NANOCOUNT } from "${subworkflowsDir}/read_count/nanocount" addParams(LABEL: 'big_mem', EXTRAPARS: progPars["counting--nanocount"], OUTPUT:outputCounts)
include { COUNT_AND_ANNO as HTSEQ_COUNT } from "${subworkflowsDir}/read_count/htseq" addParams(CONTAINER:"biocorecrg/htseq:30e9e9c", EXTRAPARS: progPars["counting--htseq"], OUTPUT:outputCounts, LABEL:'big_cpus')
include { GET_VERSION as HTSEQ_VER } from "${subworkflowsDir}/read_count/htseq" addParams(CONTAINER:"biocorecrg/htseq:30e9e9c")

include { GET_VERSION as BAMBU_VER } from "${subworkflowsDir}/assembly/bambu"
include { ASSEMBLE as BAMBU_ASSEMBLE } from "${subworkflowsDir}/assembly/bambu" addParams(EXTRAPARS: progPars["discovery--bambu"], OUTPUT:outputAssembly, LABEL:'big_mem_cpus')

include { GET_VERSION as ISOQUANT_VER } from "${subworkflowsDir}/assembly/isoquant"
include { ASSEMBLE as ISOQUANT_ASSEMBLE } from "${subworkflowsDir}/assembly/isoquant" addParams(EXTRAPARS: progPars["discovery--isoquant"], OUTPUT:outputAssembly, LABEL:'big_mem_cpus', CONTAINER:'quay.io/biocontainers/isoquant:3.2.0--hdfd78af_0')

include { REPORT as MULTIQC; GET_VERSION as MULTIQC_VER } from "${subworkflowsDir}/reporting/multiqc" addParams(EXTRAPARS: "-c ${config_report.getName()}", OUTPUT:outputMultiQC)
include { concatenateFastQFiles} from "${local_modules}" addParams(OUTPUT:outputFastq)
include { MinIONQC} from "${local_modules}" addParams(OUTPUT:outputQual, LABEL: 'big_mem_cpus')
include { bam2stats; countStats; joinCountStats; joinAlnStats} from "${local_modules}"
include { cleanFile as fastqCleanFile; cleanFile as bamCleanFile; cleanFile as fast5CleanFile} from "${local_modules}"
include { AssignReads} from "${local_modules}" addParams(OUTPUT:outputAssigned)
include { bam2Cram } from "${local_modules}" addParams(OUTPUT:outputCRAM, LABEL: 'big_cpus_ignore')
include { getFast5 } from "${local_modules}"


/*
* Wrapper for FILTERING
*/
workflow SEQFILTER {
    take:
        raw_bc_fastq

    main:
    // Optional fastq filtering
    switch(params.filtering) {
        case "nanofilt":
            bc_fastq = NANOFILT_FILTER(raw_bc_fastq)
            break;
        case "nanoq":
            bc_fastq = NANOQ_FILTER(raw_bc_fastq)
            break;
        default:
            bc_fastq = raw_bc_fastq
            break;
    }

    emit:
        out = bc_fastq

}

/*
* Wrapper for MAPPING
*/
workflow MAPPING {

    take:
        bc_fastq

    main:

    // Perform mapping on fastq files
    if (params.mapping == "NO") {
        stats_aln = Channel.empty()
        stats_counts = Channel.empty()
        sorted_alns = Channel.empty()
        nanoplot_qcs = Channel.empty()
        aln_indexes = Channel.empty()
        aln_reads = Channel.empty()
    }
    else {
        switch(params.mapping) {
            case "graphmap":
            //GRAPHMAP cannot align RNA, WE NEED TO CONVERT
             dna_bc_fastq = RNA2DNA(bc_fastq)
             aln_reads = GRAPHMAP(dna_bc_fastq, reference)
            break
            case "graphmap2":
             aln_reads = GRAPHMAP2(bc_fastq, reference)
            break
            case "minimap2":
             aln_reads = MINIMAP2(bc_fastq, reference)
            break
            case "bwa":
             aln_reads = BWA(reference, bc_fastq)
            break
            default:
            break

        }
    }

    emit:
        out = aln_reads
}

/*
* Wrapper for COUNTING
*/
workflow COUNTING {

    take:
        sorted_alns
        aln_indexes

    main:

    // OPTIONAL Perform COUNTING / ASSIGNMENT
    if (params.counting == "nanocount" && params.ref_type == "transcriptome") {
        read_counts = NANOCOUNT(sorted_alns.join(aln_indexes))
        assignments = AssignReads(sorted_alns, "nanocount")
        stat_counts = countStats(assignments)
        stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
    }
    else if (params.counting == "htseq" && params.ref_type == "genome") {
        htseq_out = HTSEQ_COUNT(annotation, sorted_alns.join(aln_indexes))
        read_counts = htseq_out.counts
        assignments = AssignReads(htseq_out.bam, "htseq")
        stat_counts = countStats(assignments)
        stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
    } else if (params.counting == "NO") {
        stats_counts = Channel.empty()
    } else {
        println "ERROR ################################################################"
        println "${params.counting} is not compatible with ${params.ref_type}"
        println "htseq requires a genome as reference and an annotation in GTF"
        println "nanocount requires a transcriptome as a reference"
        println "ERROR ################################################################"
        println "Exiting ..."
        System.exit(0)
    }


    emit:
        stats_counts = stats_counts

}


/*
* Wrapper for ASSEMBLY
*/
workflow ASSEMBLY {

    take:
        sorted_alns
        reference
        annotation

    main:

    if (params.discovery == "bambu" && params.ref_type == "genome"){
        sorted_alns.map{
            [it[1]]
        }.collect().map{
            ["assembly", it]
        }.set{data_to_bambu}
        BAMBU_ASSEMBLE(reference, annotation, data_to_bambu)
    } else if (params.discovery == "isoquant" && params.ref_type == "genome"){
        aln_indexes.map{
            [it[1]]
        }.collect().map{
            ["assembly", it]
        }.set{ixd_4_isoquant}

        sorted_alns.map{
            [it[1]]
        }.collect().map{
            ["assembly", it]
        }.join(ixd_4_isoquant).set{data_to_isoquant}

        ISOQUANT_ASSEMBLE(reference, annotation, data_to_isoquant)
    } else if (params.discovery == "NO") {
    } else {
        println "ERROR ################################################################"
        println "${params.discovery} is not compatible with ${params.ref_type}"
        println "bambu requires a genome as reference and an annotation in GTF"
        println "ERROR ################################################################"
        println "Exiting ..."
        System.exit(0)
    }
}

workflow BASECALL_MOP {

    take:
        input_fast5

    main:
	if (params.basecalling != "NO" ) {
		outbc = BASECALL(input_fast5)
		basecalled_fastq = outbc.basecalled_fastq
	} else {
		basecalled_fast5 = input_fast5
		basecalling_stats = channel.empty()
		basecalled_fastq = channel.empty()
	}

    emit:
	basecalled_fastq
	basecalling_stats
	basecalled_fast5 
	
}


workflow {

    analysis_type = checkInput(params.fast5, params.fastq)

    switch(analysis_type) {
        // INPUT IS RAW NANOPORE DATA
        case "fast5":
        fast5_4_analysis = getFast5(params.fast5)
        // BASECALL ONLY
        if (params.demultiplexing == "NO" ) {
            outbc = BASECALL(fast5_4_analysis)
            basecalled_fastq = outbc.basecalled_fastq
            bc_stats = reshapeSamples(outbc.basecalling_stats)
        }
        else { // BASECALL AND DEMULTIPLEX
            switch(params.demultiplexing) {
                case "deeplexicon":
                case "seqtagger":
                outbc = BASECALL(fast5_4_analysis)
                demux = DEMULTIPLEX(fast5_4_analysis, outbc.basecalled_fastq)
                demufq = demux.demultiplexed_fastq
                bc_stats = reshapeSamples(outbc.basecalling_stats)
                bc_demux_stats = reshapeSamples(demux.demultiplexed_tsv).groupTuple()
                break;
        
                case "guppy":
                case "readucks":
                outbc = BASECALL_DEMULTIPLEX(fast5_4_analysis)
                demufq = outbc.demultiplexed_fastqs
                bc_stats = reshapeSamples(outbc.basecalling_stats)
                bc_demux_stats = reshapeSamples(outbc.basecalling_stats).groupTuple()
                outbc.basecalled_fast5.view()
                break;
        
                case "dorado":
                break;
            }
        
        bc_stats = reshapeSamples(outbc.basecalling_stats)

        reshapedPrefiltDemufq = demufq.transpose().map{
            [it[1].name.replace(".fastq.gz", "").replace(".fq.gz", ""), it[1] ]
        }

        // FILTER BARCODES FOR FASTQ
        if (params.barcodes != "") {
            log.info "*********************************************************************"
            log.info "*************** Selecting only the requested barcodes ***************"
            log.info "*********************************************************************"
            basecalled_fastq = filterPerBarcodes(barcodes_to_include, reshapedPrefiltDemufq)
        } else {
           basecalled_fastq = reshapedPrefiltDemufq
        }

        basecalled_fastq.ifEmpty{exit 1, "NO COMBINATION SAMPLE---BARCODEID WAS FOUND\nPLEASE CHECK YOUR BARCODE LIST\nENDING NOW, BYE!!!"}
        
        // DEMULTI FAST5
        if (demulti_fast5_opt == "ON") {
            outbc.basecalled_fast5.view()
            basecalled_fast5 = reshapeSamples(outbc.basecalled_fast5).transpose().groupTuple()
            if (params.barcodes == "") {
                DEMULTI_FAST5(bc_demux_stats, basecalled_fast5)
            } else {
                // FILTER BARCODES FOR FAST5
                DEMULTI_FAST5_FILTER(bc_demux_stats, basecalled_fast5, barcodes_to_include)
            }
        } 
    }  
  
    // Perform MinIONQC on basecalling stats
    basecall_qc = MinIONQC(bc_stats.groupTuple())
    multiqc_data = multiqc_data.mix(basecall_qc.QC_folder.map{it[1]})

    // SEQUENCE FILTERING
    bc_fastq = SEQFILTER(basecalled_fastq).out

    // SEQUENCE ALIGNMENT
    alns = MAPPING(bc_fastq).out

    // Concatenate fastq and BAM files differently depending on if demultiplexed or not
    if (params.demultiplexing == "NO" ) {
        reshaped_bc_fastq = reshapeSamples(bc_fastq)
        reshaped_aln_reads = reshapeSamples(alns)
    } else {
        reshaped_bc_fastq = reshapeDemuxSamples(bc_fastq)
        reshaped_aln_reads = reshapeDemuxSamples(alns)
    }
        
    jaln_reads = SAMTOOLS_CAT(reshaped_aln_reads.groupTuple())
    fastq_files = concatenateFastQFiles(reshaped_bc_fastq.groupTuple())
    break
    
    // INPUT IS BASECALLED SEQUENCES
    case "fastq":
        fastq_files = Channel.fromFilePairs( params.fastq , size: 1, checkIfExists: true)
        jaln_reads = MAPPING(fastq_files).out
        break
    }

    // Perform SORTING and INDEXING on bam files
    sorted_alns = SAMTOOLS_SORT(jaln_reads)
    aln_indexes = SAMTOOLS_INDEX(sorted_alns)

    // Converting BAM to CRAM and
    if (params.cram_conv == "YES") {
        good_ref = checkRef(reference)
        bam2Cram(good_ref, params.subsampling_cram, sorted_alns.join(aln_indexes))
    }

    // Perform bam2stats on sorted bams
    aln_stats = bam2stats(sorted_alns)
    stats_aln = joinAlnStats(aln_stats.map{ it[1]}.collect())

    // Perform NanoPlot on sorted bams
    nanoplot_qcs = NANOPLOT_QC(sorted_alns)

    // Perform fastqc QC on fastq
    fastqc_files = FASTQC(fastq_files)
    multiqc_data = multiqc_data.mix(fastqc_files.map{it[1]})


    stats_counts = COUNTING(sorted_alns, aln_indexes).stats_counts    
    multiqc_data = multiqc_data.mix(stats_counts)

    // REVISE THIS
    //ASSEMBLY(sorted_alns, reference, params.annotation)

    // Perform MULTIQC report
    MULTIQC(multiqc_data.collect())

    //all_ver = BAMBU_VER().mix(DEMULTIPLEX_VER()).mix(NANOQ_VER()).mix(NANOFILT_VER())
    //.mix(GRAPHMAP_VER()).mix(GRAPHMAP2_VER())
    //.mix(MINIMAP2_VER()).mix(BWA_VER()).mix(FASTQC_VER())
    //.mix(SAMTOOLS_VERSION()).mix(NANOPLOT_VER()).mix(NANOCOUNT_VER()).mix(HTSEQ_VER()).mix(MULTIQC_VER())
    //.collectFile(name: 'tool_version.txt', newLine: false, storeDir:outputMultiQC)

 }


workflow.onComplete {

    def text = final_message("MoP3")
    println text
    if (params.hook != "") {
       notify_slack(text, params.hook)
    }
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
     def msg = final_message("MoP3")
        sendMail(to: params.email, subject: "MoP3 - preprocess execution", body: msg, attach: "${outputMultiQC}/multiqc_report.html")
    }
}
