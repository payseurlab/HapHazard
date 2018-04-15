
# This script will retrieve block length information

#specify the size of a list/matrix/dataframe in advance, sims x generations x chromosomes x demes
df_row <- number_of_sims * length(list_of_demes) * length(list_of_chromosomes) * length(list_of_generations)
df_col <- 13
block_data_matrix <- matrix(ncol=df_col,nrow=df_row)
column_names <- c("Sim","Chromosome","Gen","Deme","Type","Mean_Start_Born","Mean_Length","Var_Length","Length_LCI","Length_UCI","Jun_per_cM","JPC_LCI","JPC_UCI")

index <- 1 # a counter for adding entries to the dataframe
# file path should be WD/Eddard_S/Eddard_S_0/Eddard_S_0/BL_CHR1/
for(sim in 0:(number_of_sims - 1) )
{
  for(chr in list_of_chromosomes)
  {
    for(gen in list_of_generations)
    {
      block_filename <- paste(workdir,experiment_name,"/",experiment_name,"_",sim,"/",experiment_name,"_",sim,"/","BL_CHR",chr,'/0.',gen-1,'.bla',sep="")
      if( file.exists(block_filename) )
      {
        block_data <- read.csv(block_filename,header=TRUE)
        for(deme in list_of_demes)
        {
          deme_block_data <- subset(block_data,block_data$Deme==deme)
  
          chromosome_type <- floor(mean(deme_block_data$Type))
          
          deme_block_msb <- mean(deme_block_data$Start_Born)
          deme_block_mean_length <- mean(deme_block_data$Length)
          deme_block_length_var <- var(deme_block_data$Length)
          
          BL_LCI_UCI <- conf_interval(deme_block_data$Length,0.95)
          deme_block_length_UCI <- BL_LCI_UCI[2] 
          deme_block_length_LCI <- BL_LCI_UCI[1]
            
          deme_block_mean_jun_per_cM <- 0.01 / deme_block_mean_length
          deme_block_UCI_jpcM <- 0.01 / deme_block_length_UCI
          deme_block_LCI_jpcM <- 0.01 / deme_block_length_LCI
            
          block_data_matrix[index,1] <- sim
          block_data_matrix[index,2] <- chr
          block_data_matrix[index,3] <- gen
          block_data_matrix[index,4] <- deme
          block_data_matrix[index,5] <- chromosome_type
          block_data_matrix[index,6] <- deme_block_msb
          block_data_matrix[index,7] <- deme_block_mean_length
          block_data_matrix[index,8] <- deme_block_length_var
          block_data_matrix[index,9] <- deme_block_length_LCI          
          block_data_matrix[index,10] <- deme_block_length_UCI
          block_data_matrix[index,11] <- deme_block_mean_jun_per_cM
          block_data_matrix[index,12] <- deme_block_LCI_jpcM
          block_data_matrix[index,13] <- deme_block_UCI_jpcM
          
          index <- index + 1
        }
      }
    }
  }
}

block_dataframe <- data.frame(block_data_matrix)
colnames(block_dataframe) <- column_names

Block_Data <- list(name=experiment_name,data=block_dataframe)
Block_Data_Filename <- paste(experiment_name,'.Block_data.R',sep="")
save(Block_Data, file=Block_Data_Filename)

#clean up the workspace at the end
rm(column_names, block_dataframe,block_data, block_data_matrix, index)
rm(deme_block_data,block_filename,BL_LCI_UCI,chr,chromosome_type,deme)
rm(deme_block_LCI_jpcM,deme_block_length_LCI,deme_block_length_UCI,deme_block_length_var)
rm(deme_block_msb,deme_block_mean_length,deme_block_mean_jun_per_cM,deme_block_UCI_jpcM)
rm(df_row,df_col,gen,sim)


setwd(workdir)
