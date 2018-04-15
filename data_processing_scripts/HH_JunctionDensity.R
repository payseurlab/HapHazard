
# This script will retrieve block length information

#specify the size of a list/matrix/dataframe in advance, sims x generations x chromosomes x demes
df_row <- number_of_sims * length(list_of_demes) * length(list_of_chromosomes) * length(list_of_generations)
window_positions <- seq(0,max(list_of_chr_lengths),window_size)
column_names <- c("Sim","Chromosome","Gen","Deme","Type","N_Haplos",window_positions)
df_col <- length(column_names)
block_data_matrix <- matrix(ncol=df_col,nrow=df_row)

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
          if(chr==0)
          {
            deme_block_data <- subset(deme_block_data,deme_block_data$Type==2)
          }
           
          number_of_chromosomes <- length(unique(deme_block_data$Chr))
          chromosome_type <- floor(mean(deme_block_data$Type))
          sample_size <- length(unique(deme_block_data$Haplotype))
          
          block_data_matrix[index,1] <- sim
          block_data_matrix[index,2] <- chr
          block_data_matrix[index,3] <- gen
          block_data_matrix[index,4] <- deme
          block_data_matrix[index,5] <- chromosome_type
          block_data_matrix[index,6] <- sample_size
          
          position_index <- 7
          for( pos in window_positions)
          {
            if(pos==0)
            {
              window_junctions <- subset(deme_block_data, deme_block_data$Start>pos)  
            }
            else
            {
              window_junctions <- subset(deme_block_data, deme_block_data$Start>=pos)
            }
            window_junctions <- subset(window_junctions, window_junctions$Start<(pos+window_size))
            
            block_data_matrix[index,position_index] <- length(window_junctions$Start)
            position_index <- position_index + 1
          }
          
          block_data_matrix[index,7:df_col] <- block_data_matrix[index,7:df_col] / sample_size
          
          index <- index + 1
        }
      }
      else
      {
        block_data_matrix[index,1] <- sim
        block_data_matrix[index,2] <- chr
        block_data_matrix[index,3] <- gen
        block_data_matrix[index,4] <- NA
        block_data_matrix[index,5] <- NA
        block_data_matrix[index,6:df_col] <- NA
      }
    }
  }
}

junction_dataframe <- data.frame(block_data_matrix)
colnames(junction_dataframe) <- column_names

Junction_Data <- list(name=experiment_name,data=junction_dataframe)
Junction_Data_Filename <- paste(experiment_name,'.Junction_data.R',sep="")
save(Junction_Data, file=Junction_Data_Filename)

#clean up the workspace at the end
#rm(column_names, junction_dataframe,block_data, block_data_matrix, index)
#rm(df_row,df_col,gen,sim)

setwd(workdir)
