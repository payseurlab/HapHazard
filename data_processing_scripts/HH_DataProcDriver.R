
workdir <- "~/Documents/Homer/" # The folder containing zipped experiment results

setwd(workdir)

source("HH_analysis_functions.R")

experiment_name <- "Homer" 
number_of_sims <- 1000
list_of_demes <- c(0,1,2,3,5,6,7,8,9)
list_of_generations <- c(250,500,750,1000,1250,1500,1750,2000)
list_of_chromosomes <- c(0,1)
list_of_chr_lengths <- c(1.0,1.0)

init_anc_freq <- c(0,0,0,0,0,1,1,1,1,1) # the initial frequencies of ancestry "1"
SSC_mig_mat <- matrix(nrow=10, ncol=10) # the migration matrix
SSC_mig_rates <- rep(0.01,10) # a vector of migration rates, unidirectional, or one half the total exchange between two demes
window_size <- 0.01


source("HH_ExtractExperiment.R")
source("HH_GetHet.R")
source("HH_BlockLength.R")
source("HH_JunctionDensity.R")

#clean up the work space
#rm(sim,chr,list_of_chromsomes,list_of_demes,list_of_generations,number_of_sims)
#rm(list_of_chromsomes)
