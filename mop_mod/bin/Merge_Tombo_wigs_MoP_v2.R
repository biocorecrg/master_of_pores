###Script to process Tombo data: 

##Import required libraries:
library('plyr')
library('argparse')

##Argument parser:
#Create parser object
parser <- ArgumentParser()

#Create arguments:
parser$add_argument("-stats_wig", "--Statistic_wig", nargs=1, type="character", help="Path to Epinano mismatch results.")
parser$add_argument("-Cov_WT", "--Coverage_WT_bedgraph", nargs=1, type="character", help="Path to Epinano insertion results.")
parser$add_argument("-Cov_Control", "--Coverage_control_bedgraph", nargs=1, type="character", help="Path to Epinano deletion results.")
parser$add_argument("-output", "--Output_name", type="character", help="Output(s) filenames.")

#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

##Import data: 
statistic_wig <- read.delim(args$Statistic_wig, sep=" ", stringsAsFactors = FALSE)
coverage_WT <- read.delim(args$Coverage_WT_bedgraph, sep=" ", stringsAsFactors = FALSE)
coverage_control <- read.delim(args$Coverage_control_bedgraph, sep=" ", stringsAsFactors = FALSE)
output_name <- args$Output_name

if (nrow(statistic_wig)>1 && nrow(coverage_control)>1 && nrow(coverage_WT)>1){
    
  ##Parse coverage data:
  coverage_samples <- list(coverage_WT, coverage_control)
  coverage_samples_name <- c('Sample', 'IVT')
  
  for (j in 1:length(coverage_samples)){
    #Changing colnames:
    colnames(coverage_samples[[j]]) <- c(paste('Coverage', coverage_samples_name[j], sep='_'), 'Chr', 'Position', 'Include')
    
    #Adding chromosome information per position: 
    for (row in 1:nrow(coverage_samples[[j]])){
      element <- coverage_samples[[j]][[row,1]]
      
      #If it is the start of a new chromosome, we need to capture its name and starting position: 
      if (element == 'fixedStep'){
        chrom <- sub(".*=", "", coverage_samples[[j]][row,2])
        start_position <- sub(".*=", "", coverage_samples[[j]][row,3])
        coverage_samples[[j]][[row,4]] <- FALSE
        initial_position <- TRUE
        
      } else {
        if (initial_position==TRUE){
          coverage_samples[[j]][[row,2]] <- chrom
          coverage_samples[[j]][[row,3]] <- start_position
          coverage_samples[[j]][[row,4]] <- TRUE
          
          initial_position <- FALSE
          i <- as.numeric(start_position)
          
        } else {
          i <- i + 1
          
          #Update data:
          coverage_samples[[j]][[row,2]] <- chrom
          coverage_samples[[j]][[row,3]] <- i
          coverage_samples[[j]][[row,4]] <- TRUE
        }
        
      }
    }
    
    data_filtered <- subset(coverage_samples[[j]][-c(5)], Include==TRUE)
    data_filtered$Ref_Position <- paste(data_filtered$Chr, data_filtered$Position, sep = "_")
    assign(paste(coverage_samples_name[j],'filtered', sep='_'), data_filtered[c(5,1)])
     
  }
  
  #Merge both coverage tables: 
  final_coverage <- join(Sample_filtered, IVT_filtered, by='Ref_Position')[-c(4)]
  
  
  ##Parse stadistic data:
  #Changing colnames:
  colnames(statistic_wig) <- c('Position', 'statistic', 'Chr', 'Include')
  
  #Adding chromosome information per position: 
  for (row in 1:nrow(statistic_wig)){
    element <- statistic_wig[[row,1]]
    
    #If it is the start of a new chromosome, we need to capture its name and starting position: 
    if (element == 'fixedStep' || element == 'variableStep'){
      chrom <- sub(".*=", "", statistic_wig[row,2])
      statistic_wig[[row,4]] <- FALSE
  
    } else {
      
      #Update data:
      statistic_wig[[row,3]] <- chrom
      statistic_wig[[row,4]] <- TRUE
    
    }
  }
  
  stat_filtered <- subset(statistic_wig, Include==TRUE)
  stat_filtered$Ref_Position <- paste(stat_filtered$Chr, stat_filtered$Position, sep = "_")
  statistic_filtered <- stat_filtered[c(5,3,1,2)]
  
  
  ##Merging p-value and coverage data: 
  tombo_output <- join(statistic_filtered, final_coverage, by='Ref_Position')
  
  ##Mean p-value per sliding window (5-mer):
  #Extract all chromosomes:
  unique_chr <- unique(tombo_output$Chr)
  
  #Loop over all positions for every chromosome:
  sum_pvalue_kmer <- c()
  for (single_chr in unique_chr){
    data_analysis <- subset(tombo_output[c(2,3,4)], Chr==single_chr)
    
    #Looping over the rows while checking that positions are consecutive:
    for (row in 1:nrow(data_analysis)){
      if (row<3){
        sum_pvalue_kmer <- c(sum_pvalue_kmer, NA)
        
      } else {
      #Extract 5kmer positions:  
      kmer <- as.integer(c(data_analysis[(row-2),2], data_analysis[(row-1),2], data_analysis[row,2],
                data_analysis[(row+1),2], data_analysis[(row+2),2]))
      
      test_consecutive <- rle(diff(kmer))
      are_consecutive <- all(test_consecutive$values==1)
      
      #If consecutive, calculate mean:
      if (test_consecutive$lengths==4 && are_consecutive==TRUE){
        pvalues <- as.numeric(c(data_analysis[(row-2),3], data_analysis[(row-1),3], data_analysis[row,3],
                                 data_analysis[(row+1),3], data_analysis[(row+2),3]))
        
        sum_pvalue_kmer <- c(sum_pvalue_kmer, sum(pvalues))
      } else {
        sum_pvalue_kmer <- c(sum_pvalue_kmer, NA)
      }
      
      }
      
    }
    
  }
  tombo_output$statistic_kmer <- sum_pvalue_kmer
  tombo_output_filtered <- subset(tombo_output, statistic>=0)
  colnames(tombo_output_filtered) <- c('Ref_Position', 'Chr', 'Position', 'Tombo_SiteScore', 'Coverage_Sample', 'Coverage_IVT',
                                       'Tombo_KmerScore')
  
  #Output results table:
  write.table(tombo_output_filtered, file = paste(output_name, "_Tombo_Output.tsv", sep = ""), sep = "\t", row.names=FALSE)
}

