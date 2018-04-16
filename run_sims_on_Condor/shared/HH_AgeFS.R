
# The name of the experiment to be analyzed is given on the command line

# Move to the directory containing the experiment
#setwd(paste("~/Desktop/DJS2/", expName, sep=""))


print("Running AgeFS.R")

# The following parameters tell the program which files to analyze
#chromosomes <- c(1) # the number indexing the chromosome(s) by order in the genome to be analyzed
#generations <- c(499,999) # the numbers of the generations ( t - 1 actually) of the generations that need to be analyzed
#demes <- c(0:9) # the numbers of the demes in the stepping stone cline, these are used for block length stats
# each file is named and organized in the file hierarchy by simualtion number, chromosome, and generation
for(chr in chromosomes)
{
  for(gen in generations)
  {
      # We need to construct the correct filename for each of the block files
      # occasionally, odd file extensions occur despite the file containing the correct information
      # we need to check to see which extension the file has before attempting to open it
      blockFileName <- "a"
      blockFileName1 <- paste(expName,"_",jobID,"/BL_CHR",chr,"/0.",gen,".bla", sep="")
      blockFileName2 <-  paste(expName,"_",jobID,"/BL_CHR",chr,"/0.",gen,".b0", sep="")
      blockFileName3 <-  paste(expName,"_",jobID,"/BL_CHR",chr,"/0.",gen,".bl0", sep="")
      # check to see if each file exists, print the name of the one that does,
      # then for the file that exists, read in the data, and assign that filename
      # to the blockFileName variable
      if(file.exists(blockFileName1))
      {
        print(blockFileName1)
        blockLengths <- read.csv(blockFileName1)
        blockFileName <- blockFileName1
      }
      else if (file.exists(blockFileName2))
      {
        print(blockFileName2)
        blockLengths <- read.csv(blockFileName2)
        blockFileName <- blockFileName2
      }
      else if (file.exists(blockFileName3))
      {
        print(blockFileName3)
        blockLengths <- read.csv(blockFileName3)
        blockFileName <- blockFileName3
      }
      else
      {    
        print(paste("Could not open file: ", blockFileName1, sep=""))
        blockFileName <- blockFileName1
      }
      
    Age.by.deme <- list()
    block.data <- read.csv(blockFileName)
    block.data <- subset(block.data, block.data[,9] != 0)
    Age.by.deme[[1]] <- hist(block.data[,9],breaks=gen+1, plot=FALSE)
    
    for(d in demes)
    {
      block.deme <- subset(block.data, block.data[,1] == d)
      Age.by.deme[[d+2]] <- hist(block.deme[,9], breaks=gen+1, plot=FALSE)
    }
    
    Freq.Spec <- list()
    junlist <- subset(block.data, block.data[,7] > 0)
    junlist <- subset(block.data, block.data[,7] < 1)
    junlist <- junlist[,8]
    
    junlist <- table(junlist)
    junlist <- data.frame(junlist)
    Freq.Spec[[1]] <- junlist
    
    for(d in demes)
    {
      junlist <- subset(block.data, block.data[,1] == d)
      junlist <- subset(junlist, junlist[,7] > 0)
      junlist <- subset(junlist, junlist[,7] < 1)
      junlist <- junlist[,8]
      
      junlist <- table(junlist)
      junlist <- data.frame(junlist)
      Freq.Spec[[d+2]] <- data.frame(junlist)    
    }
    
    # Save all the matrices with the data to files
    save(Age.by.deme, file=paste(expName,"_",jobID,"/",expName, "_" , "0_", chr, "_", gen, "_age.R", sep="") )
    save(Freq.Spec, file=paste(expName,"_",jobID,"/",expName, "_" , "0_", chr, "_", gen, "_frs.R", sep="") )
    
  }
}
