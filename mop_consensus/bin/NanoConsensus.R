###Script to merge results from NanoMod and generate a consensus putative modified positions list:

##Import libraries:
suppressMessages(library('plyr'))
suppressMessages(library('dplyr'))
suppressMessages(library('VennDiagram'))
suppressMessages(library('argparse'))
suppressMessages(library('ggplot2'))
suppressMessages(library('GenomicRanges'))
suppressMessages(library('stringr'))
suppressMessages(library('scales'))
suppressMessages(library('ggnewscale'))
suppressMessages(library('ggrepel'))
suppressMessages(library('gtable'))
suppressMessages(library('xfun'))


##Import accessory functions:
source('./scripts/Accessory_functions_consensusNanoMod.R')

##Argument parser:
#Create parser object
parser <- ArgumentParser()

#Define desired outputs:
#GLOBAL FEATURES:
parser$add_argument("-output", "--Output_name", type="character", help="Output(s) filenames.")
parser$add_argument("-fasta", "--Fasta_file", type="character", help="Genome fasta file.")
parser$add_argument("-ini_pos", "--Initial_position", type="integer", default=50, help="Initial position [default %(default)].")
parser$add_argument("-fin_pos", "--Final_position", type="integer", help="Final position.")
#parser$add_argument("-plot", "--Plotting", action="store_true", help="Plot significant positions for all methods.")
parser$add_argument("-chr", "--Chr", type="character", help="Character to study.")
parser$add_argument("--MZS_thr", default=5, type="double", 
                    help="Modified Z-Score threshold for all results [default %(default)]")
parser$add_argument("--NC_thr", default=5, type="double", 
                    help="NanoConsensus score threshold for all results [default %(default)]")
parser$add_argument("-exclude", "--Exclude", nargs='+', type="integer", help="Exclude these positions from the analysis (SNPs) - it will exclude the 17-mer.")
parser$add_argument("--model_score", default="global", type="character", 
                    help="Model used to calculate NanoConsensus score [default %(default)]")
parser$add_argument("--coverage", default=1, type="integer", 
                    help="Minimum coverage per position to be included in the analysis [default %(default)]")
parser$add_argument("--nanocomp_stat", default="GMM_logit_pvalue_context_2", type="character", 
                    help="Stat from Nanocompore output to be used [default %(default)]")

#EPINANO:
parser$add_argument("-Epi_Sample", "--Epinano_Sample", nargs=1, type="character", help="Path to Epinano features sample results.")
parser$add_argument("-Epi_IVT", "--Epinano_IVT", nargs=1, type="character", help="Path to Epinano features IVT results.")

#NANOPOLISH:
parser$add_argument("-NP_Sample", "--Nanopolish_Sample", nargs=1, type="character", help="Path to Nanopolish mean per position sample results.")
parser$add_argument("-NP_IVT", "--Nanopolish_IVT", nargs=1, type="character", help="Path to Nanopolish mean per position IVT results.")

#TOMBO:
parser$add_argument("-Tombo", "--Tombo_Sample", nargs=1, type="character", help="Path to Tombo pairwise comparison results.")

#NANOCOMPORE:
parser$add_argument("-Nanocomp", "--Nanocomp_Sample", nargs=1, type="character", help="Path to Nanocompore pairwise comparison results.")


#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

##Create and update log file:
write('NanoConsensus - v 1.0', file = paste("NanoConsensus_", args$Output_name,".log", sep=""))
write(paste('Analysing sample: ',args$Output_name, sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
 write(paste('Minimum coverage: ',args$coverage, sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write(paste('Transcript analysed: ',args$Chr, " - from ", args$Initial_position, " to ", args$Final_position, sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write(paste('Nanocompore stat used: ',args$nanocomp_stat, sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write(paste('Z score threshold: ',args$MZS_thr, sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write(paste('NanoConsensus score threshold: ',args$NC_thr, "*median(NanoConsensus Score)", sep=""), file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
write('Step 1: Processing data from individual softwares', file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)

##EPINANO processing: 
epinano_data <- epinano_processing(args$Epinano_Sample, args$Epinano_IVT, args$Initial_position, args$Final_position, args$MZS_thr, args$Chr, args$Exclude, args$coverage)

##NANOPOLISH processing: 
nanopolish_data <- nanopolish_processing(args$Nanopolish_Sample, args$Nanopolish_IVT, args$Initial_position, args$Final_position, args$MZS_thr, args$Chr, args$Exclude, args$coverage)

##TOMBO processing: 
tombo_data <- tombo_processing(args$Tombo_Sample, args$thr_tombo_pos, args$thr_tombo_kmer, args$Initial_position, args$Final_position, args$MZS_thr, args$Chr, args$Exclude, args$coverage)

##NANOCOMPORE processing:
nanocompore_data <- nanocomp_processing(args$Nanocomp_Sample, args$nanocomp_metric, args$thr_nanocomp, args$Initial_position, args$Final_position, args$MZS_thr, args$Chr, args$Exclude, args$nanocomp_stat)

##DATA PROCESSING:
#Generate list with all positions and significant positions respectively:
list_plotting <- list(epinano_data[[1]], nanopolish_data[[1]], tombo_data[[1]], nanocompore_data[[1]])
list_significant <- list(epinano_data[[2]], nanopolish_data[[2]], tombo_data[[2]], nanocompore_data[[2]])

#Create Z-Scores plotting object:
write('Step 2: Plotting ZScores from individual softwares', file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
barplot_4soft <- barplot_plotting(list_plotting, list_significant, args$Output_name, args$MZS_thr, args$Autoscaling, args$Initial_position, args$Final_position)

##Analysis of SIGNIFICANT POSITIONS across methods:
write('Step 3: Overlapping analysis and generation of Venn diagram', file = paste("NanoConsensus_", args$Output_name,".log", sep=""), append = T)
analysis_significant_positions(list_significant, list_plotting, args$Fasta_file, args$Output_name,  args$Initial_position, args$Final_position, args$MZS_thr, args$NC_thr, args$model_score, barplot_4soft)

