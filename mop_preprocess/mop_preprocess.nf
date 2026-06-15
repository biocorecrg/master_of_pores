#!/usr/bin/env nextflow

nextflow.enable.dsl=2



/*
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '4.1'

params.help            = false
params.resume          = false

def local_modules = file("${projectDir}/../local_modules.nf")

include { colorCodes } from local_modules

def colors = colorCodes()

log.info """

${colors.yellow}${colors.bold}====================================================
╔╦╗╔═╗╔═╗  ╔═╗┬─┐┌─┐┌─┐┬─┐┌─┐┌─┐┌─┐┌─┐┌─┐
║║║║ ║╠═╝  ╠═╝├┬┘├┤ ├─┘├┬┘│ ││  ├┤ └─┐└─┐
╩ ╩╚═╝╩    ╩  ┴└─└─┘┴  ┴└─└─┘└─┘└─┘└─┘└─┘
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣷⣶⣤⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⣿⣿⣿⣿⣷⡒⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⣿⣿⣿⣿⣿⣆⠙⡄⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣤⣤⣤⣤⣤⣤⣤⣤⣤⠤⢄⡀⠀⠀⣿⣿⣿⣿⣿⣿⡆⠘⡄⠀⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⢿⣿⣿⣿⣿⣿⣿⣿⣦⡈⠒⢄⢸⣿⣿⣿⣿⣿⣿⡀⠱⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣿⣿⣿⣿⣿⣿⣿⣦⠀⠱⣿⣿⣿⣿⣿⣿⣇⠀⢃⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢿⣿⣿⣿⣿⣿⣿⣷⡄⣹⣿⣿⣿⣿⣿⣿⣶⣾⣿⣶⣤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣴⣶⣿⣭⣍⡉⠙⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⢀⣠⣶⣿⣿⣿⣿⣿⣿⣿⣿⣷⣦⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀⠀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠉⠉⠛⠻⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡷⢂⣓⣶⣶⣶⣶⣤⣤⣄⣀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⣿⣿⠟⢀⣴⢿⣿⣿⣿⠟⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠛⠋⠉⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠤⠤⠤⠤⠙⣻⣿⣿⣿⣿⣿⣿⣾⣿⣿⡏⣠⠟⡉⣾⣿⣿⠋⡠⠊⣿⡟⣹⣿⢿⣿⣿⣿⠿⠛⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣤⣶⣤⣭⣤⣼⣿⢛⣿⣿⣿⣿⣻⣿⣿⠇⠐⢀⣿⣿⡷⠋⠀⢠⣿⣺⣿⣿⢺⣿⣋⣉⣉⣩⣴⣶⣤⣤⣄⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠛⠻⠿⣿⣿⣿⣇⢻⣿⣿⡿⠿⣿⣯⡀⠀⢸⣿⠋⢀⣠⣶⠿⠿⢿⡿⠈⣾⣿⣿⣿⣿⡿⠿⠛⠋⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠻⢧⡸⣿⣿⣿⠀⠃⠻⠟⢦⢾⢣⠶⠿⠏⠀⠰⠀⣼⡇⣸⣿⣿⠟⠉⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣴⣾⣶⣽⣿⡟⠓⠒⠀⠀⡀⠀⠠⠤⠬⠉⠁⣰⣥⣾⣿⣿⣶⣶⣷⡶⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠉⠹⠟⣿⣿⡄⠀⠀⠠⡇⠀⠀⠀⠀⠀⢠⡟⠛⠛⠋⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⠋⠹⣷⣄⠀⠐⣊⣀⠀⠀⢀⡴⠁⠣⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣤⣀⠤⠊⢁⡸⠀⣆⠹⣿⣧⣀⠀⠀⡠⠖⡑⠁⠀⠀⠀⠑⢄⣀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⣦⣶⣿⣿⣟⣁⣤⣾⠟⠁⢀⣿⣆⠹⡆⠻⣿⠉⢀⠜⡰⠀⠀⠈⠑⢦⡀⠈⢾⠑⡾⠲⣄⠀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⣀⣤⣶⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠖⠒⠚⠛⠛⠢⠽⢄⣘⣤⡎⠠⠿⠂⠀⠠⠴⠶⢉⡭⠃⢸⠃⠀⣿⣿⣿⠡⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⡤⠶⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣋⠁⠀⠀⠀⠀⠀⢹⡇⠀⠀⠀⠀⠒⠢⣤⠔⠁⠀⢀⡏⠀⠀⢸⣿⣿⠀⢻⡟⠑⠢⢄⡀⠀⠀⠀⠀
⠀⠀⠀⠀⢸⠀⠀⠀⡀⠉⠛⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣄⣀⣀⡀⠀⢸⣷⡀⣀⣀⡠⠔⠊⠀⠀⢀⣠⡞⠀⠀⠀⢸⣿⡿⠀⠘⠀⠀⠀⠀⠈⠑⢤⠀⠀
⠀⠀⢀⣴⣿⡀⠀⠀⡇⠀⠀⠀⠈⣿⣿⣿⣿⣿⣿⣿⣿⣝⡛⠿⢿⣷⣦⣄⡀⠈⠉⠉⠁⠀⠀⠀⢀⣠⣴⣾⣿⡿⠁⠀⠀⠀⢸⡿⠁⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀
⠀⢀⣾⣿⣿⡇⠀⢰⣷⠀⢀⠀⠀⢹⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣦⣭⣍⣉⣉⠀⢀⣀⣤⣶⣾⣿⣿⣿⢿⠿⠁⠀⠀⠀⠀⠘⠀⠀⠀⠀⠀⠀⠀⠀⠀⡰⠉⢦⠀
⢀⣼⣿⣿⡿⢱⠀⢸⣿⡀⢸⣧⡀⠀⢿⣿⣿⠿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡭⠖⠁⠀⡠⠂⠀⠀⠀⠀⠀⠀⠀⠀⢠⠀⠀⠀⢠⠃⠀⠈⣀
⢸⣿⣿⣿⡇⠀⢧⢸⣿⣇⢸⣿⣷⡀⠈⣿⣿⣇⠈⠛⢿⣿⣿⣿⣿⣿⣿⠿⠿⠿⠿⠿⠿⠟⡻⠟⠉⠀⠀⡠⠊⠀⢠⠀⠀⠀⠀⠀⠀⠀⠀⣾⡄⠀⢠⣿⠔⠁⠀⢸
⠈⣿⣿⣿⣷⡀⠀⢻⣿⣿⡜⣿⣿⣷⡀⠈⢿⣿⡄⠀⠀⠈⠛⠿⣿⣿⣿⣷⣶⣶⣶⡶⠖⠉⠀⣀⣤⡶⠋⠀⣠⣶⡏⠀⠀⠀⠀⠀⠀⠀⢰⣿⣧⣶⣿⣿⠖⡠⠖⠁
⠀⣿⣿⣷⣌⡛⠶⣼⣿⣿⣷⣿⣿⣿⣿⡄⠈⢻⣷⠀⣄⡀⠀⠀⠀⠈⠉⠛⠛⠛⠁⣀⣤⣶⣾⠟⠋⠀⣠⣾⣿⡟⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣿⠷⠊⠀⢰⠀
⢰⣿⣿⠀⠈⢉⡶⢿⣿⣿⣿⣿⣿⣿⣿⣿⣆⠀⠙⢇⠈⢿⣶⣦⣤⣀⣀⣠⣤⣶⣿⣿⡿⠛⠁⢀⣤⣾⣿⣿⡿⠁⠀⠀⠀⠀⠀⠀⠀⣸⣿⡿⠿⠋⠙⠒⠄⠀⠉⡄
⣿⣿⡏⠀⠀⠁⠀⠀⠀⠉⠉⠙⢻⣿⣿⣿⣿⣷⡀⠀⠀⠀⠻⣿⣿⣿⣿⣿⠿⠿⠛⠁⠀⣀⣴⣿⣿⣿⣿⠟⠀⠀⠀⠀⠀⠀⠀⠀⢠⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠰
====================================================
BIOCORE@CRG Master of Pores 4. Preprocessing - N F  ~  version ${version}
====================================================
${colors.reset}

${colors.bold}Input${colors.reset}
----------------------------------------------------${colors.reset}
${colors.green}pod5${colors.reset}                      : ${params.pod5}
${colors.green}fastq${colors.reset}                     : ${params.fastq}

${colors.bold}Reference${colors.reset}
----------------------------------------------------${colors.reset}
${colors.green}reference${colors.reset}                  : ${params.reference}
${colors.green}annotation${colors.reset}                 : ${params.annotation}
${colors.green}ref_type${colors.reset}                   : ${params.ref_type}

${colors.bold}Output${colors.reset}
----------------------------------------------------${colors.reset}
${colors.green}output${colors.reset}                    : ${params.output}
${colors.green}email${colors.reset}                     : ${params.email}
${colors.green}slackhook${colors.reset}                 : ${params.slackhook}

${colors.bold}Actions
----------------------------------------------------${colors.reset}
${colors.green}basecalling${colors.reset}               : ${params.basecalling}
${colors.green}demultiplexing${colors.reset}            : ${params.demultiplexing}
${colors.green}demulti_pod5${colors.reset}              : ${params.demulti_pod5}
${colors.green}filtering${colors.reset}                 : ${params.filtering}
${colors.green}mapping${colors.reset}                   : ${params.mapping}
${colors.green}counting${colors.reset}                  : ${params.counting}
${colors.green}discovery${colors.reset}                 : ${params.discovery}
${colors.green}cram_conv${colors.reset}                 : ${params.cram_conv}
${colors.green}subsampling_cram${colors.reset}          : ${params.subsampling_cram}

${colors.bold}Advanced${colors.reset}
----------------------------------------------------${colors.reset}
${colors.green}granularity${colors.reset}               : ${params.granularity}
${colors.green}barcodes${colors.reset}                  : ${params.barcodes}
${colors.green}GPU${colors.reset}                       : ${params.GPU}

${colors.bold}====================================================${colors.reset}

"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// include functions, outdirs from other files
evaluate(new File("../outdirs.nf"))
//def local_modulesDir = "${projectDir}/../local"
def subworkflowsDir = "${projectDir}/../BioNextflow/subworkflows"
def workflowsDir = "${projectDir}/../BioNextflow/workflows"
joinScript = file("${projectDir}/bin/join.r")

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



// Sanitize parameters
if (params.mapping != "NO") {
    reference = file(params.reference)
    if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
	cram_conv = params.cram_conv
	counting = params.counting
	discovery = params.discovery
} else {
    cram_conv = "NO"
    reference = ""
    counting = "NO"
    discovery = "NO"
}

if (params.demultiplexing == "NO" ) {
	demulti_pod5_opt = "NO"
} else {
	demulti_pod5_opt = params.demulti_pod5
}

// Get models
demux_models = ""

switch(params.demultiplexing) {
    case "seqtagger":
        demux_models = "${projectDir}/seqtagger_models/"
    break;
    case "seqtagger-trna":
        demux_models = "${projectDir}/seqtagger_tRNA_models/"
    break;
    case "dorado":
        demux_models = "${projectDir}/seqtagger_models/"
    break;
}

dorado_models = "${projectDir}/dorado_models/"

// check GPU usage.
if (params.GPU != "LOCAL" && params.GPU != "cuda11" && params.GPU != "cuda10" && params.GPU != "OFF" && params.GPU != "ON") exit 1, "Please specify cuda11, cuda10, ON, LOCAL or OFF if GPU processors are available. ON is legacy for cuda10"

def gpu = (params.GPU != 'OFF') ? 'ON' : 'OFF'
gpu_bc = (params.GPU == 'LOCAL') ? 'LOCAL' : gpu


// CHECK INCOMPATIBILITIES AMONG PARAMETERS
annotation = ""
if (params.ref_type == "genome") {
    if (params.annotation) {
        annotation = file(params.annotation)
        if( !annotation.exists() ) exit 1, "Missing annotation file: ${params.annotation}!"
    }
}

include {final_message; notify_slack } from "${subworkflowsDir}/global_functions.nf"

include { checkInput; filterPerBarcodes; get_barcode_list; RNA2DNA; parseFinalSummary; checkTools; reshapeSamples; reshapeDemuxSamples; checkRef; getParameters; homogenizeVals } from "${local_modules}"

def demulti_pod5_opt = homogenizeVals(demulti_pod5_opt)
def basecall_label = (params.GPU != 'OFF' ? 'basecall_gpus' : 'big_cpus')

if (demulti_pod5_opt == "ON" && params.pod5) {
	log.info """${colors.green}DEMULTIPLEXING POD5${colors.reset}
	"""
}


//def output_bc = (demulti_pod5_opt == 'ON' ? '' : outputFast5)
def output_bc = ''
//def outputMinionQC = (demulti_pod5_opt == 'ON' ? '': outputQual)

def basecalling = params.basecalling


// GET PROGRAM PARS AND VERIFY
def tools = [:]
tools["basecalling"] = homogenizeVals(basecalling)
tools["demultiplexing"] = homogenizeVals(params.demultiplexing)
tools["mapping"] = homogenizeVals(params.mapping)
tools["filtering"] = homogenizeVals(params.filtering)
tools["counting"] = homogenizeVals(params.counting)
tools["discovery"] = homogenizeVals(params.discovery)

// Remove basecalling and demultiplexing in case of fastq input
if (!params.fastq) params.fastq = ""

if(params.pod5 == "" && params.fastq != "") {
    basecalling = "NO"
    tools["basecalling"] = "NO"
    tools["demultiplexing"] = "NO"
}

progPars = params.progPars
checkTools(tools, params.progPars)


// Create a channel for excluded ids
//barcodes_to_include = get_barcode_list(params.barcodes)

progPars["basecalling"]["dorado-mod"] = progPars["basecalling"]["dorado-mod"] + " --emit-moves"

def basecaller_pars = ["dorado" : progPars["basecalling"]["dorado"],  "dorado-duplex" : progPars["basecalling"]["dorado"], "dorado-mod" : progPars["basecalling"]["dorado-mod"] ]
def demux_pars = [ "dorado" : progPars["basecalling"]["dorado"] + " " + progPars["demultiplexing"]["dorado"], "seqtagger":  progPars["demultiplexing"]["seqtagger"], "seqtagger-trna":  progPars["demultiplexing"]["seqtagger-trna"] ]
def mapping_pars = ["bwa": progPars["mapping"]["bwa"], "winnowmap": progPars["mapping"]["winnowmap"] + " -y",
				"graphmap2": progPars["mapping"]["graphmap2"], "minimap2": progPars["mapping"]["minimap2"] + " -y --MD",
				"graphmap": progPars["mapping"]["graphmap"]
				]

def dem_cont = ""
if (params.demultiplexing == "seqtagger-trna") {
	dem_cont = "lpryszcz/seqtagger:1.1a"
}

// INCLUDE WORKFLOWS
include { BASECALL } from "${workflowsDir}/basecaller" addParams(gpu: gpu_bc, output: output_bc, label: basecall_label, label2:'big_cpus', type:basecalling ,  extrapars: basecaller_pars[basecalling], models: dorado_models )
include { DEMULTIPLEX } from "${workflowsDir}/demultiplexer.nf" addParams(container: dem_cont, gpu: gpu, output: output_bc, label: basecall_label, type:params.demultiplexing , extrapars: demux_pars[params.demultiplexing], models: demux_models )
include { BASECALL_DEMULTIPLEX } from "${workflowsDir}/basecaller_demultiplexer.nf" addParams(gpu: gpu_bc, output: output_bc, label: basecall_label, label2:'big_cpus', type:params.demultiplexing , extrapars: demux_pars[params.demultiplexing], models: dorado_models  )

//big_cpus_retry
include { DEMULTI_POD5; DEMULTI_POD5_FILTER } from "${subworkflowsDir}/misc/demulti_pod5" addParams(OUTPUT: outputFast5, OUTPUTST: outputQual, LABEL: 'big_cpus_retry', TYPE: params.demultiplexing)
include { ALIGN } from "${workflowsDir}/aligner.nf" addParams(output: output_bc, label: 'big_mem_cpus', type: params.mapping , extrapars: mapping_pars[params.mapping] )


// INCLUDE MODULES
include { GET_VERSION as SEQKIT_VER; FILTER as SEQKIT_FILTER} from "${subworkflowsDir}/trimming/seqkit" addParams(EXTRAPARS: progPars["filtering"]["seqkit"])
include { GET_VERSION as NANOQ_VER; FILTER as NANOQ_FILTER} from "${subworkflowsDir}/trimming/nanoq" addParams(EXTRAPARS: progPars["filtering"]["nanoq"])
include { REPORT as NANOQ_REPORT} from "${subworkflowsDir}/trimming/nanoq" addParams(EXTRAPARS: "-t 5 -vvv")

include { SORT as SAMTOOLS_SORT } from "${subworkflowsDir}/misc/samtools" addParams(LABEL: 'big_cpus_retry', OUTPUT:outputMapping)
include { INDEX as SAMTOOLS_INDEX } from "${subworkflowsDir}/misc/samtools" addParams(OUTPUT:outputMapping)
include { GET_VERSION as SAMTOOLS_VERSION; CAT as SAMTOOLS_CAT } from "${subworkflowsDir}/misc/samtools"
include { QC as NANOSTAT_QC } from "${subworkflowsDir}/qc/nanostat" addParams(LABEL: 'big_cpus_ignore')
include { GET_VERSION as NANOCOUNT_VER } from "${subworkflowsDir}/read_count/nanocount"
include { COUNT as NANOCOUNT } from "${subworkflowsDir}/read_count/nanocount" addParams(LABEL: 'big_mem', EXTRAPARS: progPars["counting"]["nanocount"], OUTPUT:outputCounts)
include { COUNT_AND_ANNO as HTSEQ_COUNT } from "${subworkflowsDir}/read_count/htseq" addParams(CONTAINER:"biocorecrg/htseq:30e9e9c", EXTRAPARS: progPars["counting"]["htseq"], OUTPUT:outputCounts, LABEL:'big_cpus')
include { GET_VERSION as HTSEQ_VER } from "${subworkflowsDir}/read_count/htseq" addParams(CONTAINER:"biocorecrg/htseq:30e9e9c")

include { GET_VERSION as BAMBU_VER } from "${subworkflowsDir}/assembly/bambu"
include { ASSEMBLE as BAMBU_ASSEMBLE } from "${subworkflowsDir}/assembly/bambu" addParams(EXTRAPARS: progPars["discovery"]["bambu"], OUTPUT:outputAssembly, LABEL:'big_mem_cpus')

include { GET_VERSION as ISOQUANT_VER } from "${subworkflowsDir}/assembly/isoquant"
include { ASSEMBLE as ISOQUANT_ASSEMBLE } from "${subworkflowsDir}/assembly/isoquant" addParams(EXTRAPARS: "--data_type nanopore  " + progPars["discovery"]["isoquant"], OUTPUT:outputAssembly, LABEL:'big_time_cpus', CONTAINER:'quay.io/biocontainers/isoquant:3.2.0--hdfd78af_0')

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
        case "seqfilter":
            bc_fastq = SEQKIT_FILTER(raw_bc_fastq)
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

workflow MAPPING_MOP {

    take:
        bc_fastq
        reference

    main:

    // Perform mapping on fastq files
    if (params.mapping == "NO") {
        aln_reads = Channel.value()
    }
    else {
        if(params.mapping == "graphmap") {
            //GRAPHMAP cannot align RNA, WE NEED TO CONVERT
             bc_fastq = RNA2DNA(bc_fastq)
         }
	aln_reads = ALIGN(bc_fastq, reference).out
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
        aln_indexes
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
	if (basecalling != "NO" ) {
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

    /* DEFAULT EMPTY VALS
    stats_counts = Channel.value()
    sorted_alns = Channel.value()
    nanoplot_qcs = Channel.value()
    aln_indexes = Channel.value()
    alns = Channel.value()
    aln_stats = Channel.value()
    */

    analysis_type = checkInput(params.pod5, params.fastq)

    switch(analysis_type) {
        // INPUT IS RAW NANOPORE DATA
        case "fast5":
        pod5_4_analysis = getFast5(params.pod5)
        // BASECALL ONLY
        if (params.demultiplexing == "NO" ) {
            outbc = BASECALL(pod5_4_analysis)
            basecalled_fastq = outbc.basecalled_fastq
            bc_stats = reshapeSamples(outbc.basecalling_stats)
        }
        else { // BASECALL AND DEMULTIPLEX

            switch(params.demultiplexing) {
                case "seqtagger":
                case "seqtagger-trna":
                	outbc = BASECALL(pod5_4_analysis)
                	demux = DEMULTIPLEX(pod5_4_analysis, outbc.basecalled_fastq)
                	demufq = demux.demultiplexed_fastq
                	bc_stats = reshapeSamples(outbc.basecalling_stats)
                	bc_demux_stats = reshapeSamples(demux.demultiplexed_tsv).groupTuple()
                	break;
                case "dorado":
                	outbc = BASECALL_DEMULTIPLEX(pod5_4_analysis)
                	demufq = outbc.demultiplexed_fastqs
                	bc_demux_stats = reshapeSamples(outbc.basecalling_stats).groupTuple()
                	break;
                default:
			        println "ERROR ################################################################"
        			println "${params.demultiplexing} is not supported!!!"
        			println "ERROR ################################################################"
        			println "Exiting ..."
        			System.exit(0)
                	break;
            }
		bc_stats = reshapeSamples(outbc.basecalling_stats)
		reshapedPrefiltDemufq = demufq.transpose().map{
			[it[1].name.replace(".fastq.gz", "").replace(".fq.gz", ""), it[1] ]
		}

		// FILTER BARCODES FOR FASTQ
		if (params.barcodes) {
                barcodes_to_include = get_barcode_list(params.barcodes)
				log.info "*********************************************************************"
				log.info "*************** Selecting only the requested barcodes ***************"
				log.info "*********************************************************************"
				basecalled_fastq = filterPerBarcodes(barcodes_to_include, reshapedPrefiltDemufq)
			} else {
			   basecalled_fastq = reshapedPrefiltDemufq
			}

			basecalled_fastq.ifEmpty{exit 1, "NO COMBINATION SAMPLE---BARCODEID WAS FOUND\nPLEASE CHECK YOUR BARCODE LIST\nENDING NOW, BYE!!!"}

			// DEMULTI POD5. POD5 are not basecalled so they just need to be split
			if (demulti_pod5_opt == "ON") {
				grouped_pod5 = reshapeSamples(pod5_4_analysis)
				if (!params.barcodes) {
					DEMULTI_POD5(bc_demux_stats, grouped_pod5)
				} else {
					// FILTER BARCODES FOR POD5
                	my_barcodes = get_barcode_list(params.barcodes).map{
                		id = it.split("---")[0]
                		[id, it]
                	}.collectFile { item ->
        				[ "${item[0]}___barcodes.txt", item[1] + '\n' ]
    				}.map {
    					id = "${it.getSimpleName()}".split("___")[0]
    					[id, it]
    				}

					DEMULTI_POD5_FILTER(bc_demux_stats, grouped_pod5, my_barcodes)
				}
			}
		}

		// SEQUENCE FILTERING
		bc_fastq = SEQFILTER(basecalled_fastq).out
		alns = MAPPING_MOP(bc_fastq, reference).out

		// Concatenate fastq and BAM files differently depending on if demultiplexed or not
		if (params.demultiplexing == "NO" ) {
			reshaped_bc_fastq = reshapeSamples(bc_fastq)
	    	if (params.mapping != "NO") {
				reshaped_aln_reads = reshapeSamples(alns)
			}
		} else {
			reshaped_bc_fastq = reshapeDemuxSamples(bc_fastq)
	    	if (params.mapping != "NO") {
				reshaped_aln_reads = reshapeDemuxSamples(alns)
			}
		}

    	if (params.mapping != "NO") {
    		jaln_reads = SAMTOOLS_CAT(reshaped_aln_reads.groupTuple())
  	  		//aln_indexes = SAMTOOLS_INDEX(jaln_reads)
    	}
		fastq_files = concatenateFastQFiles(reshaped_bc_fastq.groupTuple())
		break;

		// INPUT IS BASECALLED SEQUENCES
		case "fastq":
			fastq_files = Channel.fromFilePairs( params.fastq , size: 1, checkIfExists: true)
			jaln_reads = MAPPING_MOP(fastq_files, reference).out
			break;
		}

    // Perform SORTING and INDEXING on bam files
    if (params.mapping != "NO") {
  	  sorted_alns = SAMTOOLS_SORT(jaln_reads)
  	  aln_indexes = SAMTOOLS_INDEX(sorted_alns)
    // Perform bam2stats on sorted bams
      aln_stats = bam2stats(sorted_alns)
      stats_aln = joinAlnStats(aln_stats.map{ it[1]}.collect())
      multiqc_data = multiqc_data.mix(stats_aln)
    // Perform NanoPlot on sorted bams
       nanoplot_qcs = NANOSTAT_QC(sorted_alns)
       multiqc_data = multiqc_data.mix(nanoplot_qcs.map{ it[1]})

    }

    // Converting BAM to CRAM and
    if (cram_conv == "YES") {
        good_ref = checkRef(reference)
        bam2Cram(good_ref, params.subsampling_cram, sorted_alns.join(aln_indexes))
    }

    // Perform fastqc QC on fastq
    nqreport_files = NANOQ_REPORT(fastq_files)
    multiqc_data = multiqc_data.mix(nqreport_files.map{it[1]})

	if (counting != "NO") {
        stats_counts = COUNTING(sorted_alns, aln_indexes).stats_counts
        multiqc_data = multiqc_data.mix(stats_counts)
    }

    // REVISE THIS
    if (discovery != "NO") {
	    ASSEMBLY(sorted_alns, aln_indexes, reference, params.annotation)
	}

    // Perform MULTIQC report
    MULTIQC(multiqc_data.collect())

    //all_ver = BAMBU_VER().mix(DEMULTIPLEX_VER()).mix(NANOQ_VER()).mix(NANOFILT_VER())
    //.mix(GRAPHMAP_VER()).mix(GRAPHMAP2_VER())
    //.mix(MINIMAP2_VER()).mix(BWA_VER()).mix(FALCOQC_VER())
    //.mix(SAMTOOLS_VERSION()).mix(NANOPLOT_VER()).mix(NANOCOUNT_VER()).mix(HTSEQ_VER()).mix(MULTIQC_VER())
    //.collectFile(name: 'tool_version.txt', newLine: false, storeDir:outputMultiQC)

 }




workflow.onComplete {

    def text = final_message("MOP4")
    println text
    if (params.slackhook) {
       notify_slack(text, params.slackhook)
    }
}

/*
* Mail notification
*/

if (params.email == "yourmail@yourdomain" || !params.email) {
    log.info 'Skipping the email\n'
}
else {
    log.info "Sending the email to ${params.email}\n"

    workflow.onComplete {
     def msg = final_message("MOP4")
        sendMail(to: params.email, subject: "MOP4 - preprocess execution", body: msg, attach: "${outputMultiQC}/MOP4-pipeline_multiqc_report.html")
    }
}
