
# The name of the experiment to be analyzed is given on the command line
#args <- commandArgs(TRUE)

#for(i in 1:length(args))
#{
#	eval(parse(text=args[[i]]))
#}

print("Running HH_BJ_1s.R")

# The following parameters tell the program which files to analyze
simNum <- 0 # the number associated with each run for that experiment
#chromosomes <- c(0,1,2) # the number indexing the chromosome(s) by order in the genome to be analyzed
deme <- c(0:9) # the numbers of the demes in the stepping stone cline, these are used for block length stats
chromosomes <- seq(0,length(chrLength)-2,1)
# get the lengths of each of the chromosomes that were run in the experiement
#chrLength <- c(1,1,0)


#winSize <- 0.01 # the length in Morgans of the windows used to count junctions
# the demes for which genomic junction densities will be calculated
sampSize <- 10

# each file is named and organized in the file hierarchy by simualtion number, chromosome, and generation
for(chr in chromosomes)
{
  for(gen in generations)
  {
    # initialize matrices to store block length and junction data
    # block length is stored with a row for each simulation, and a column for each deme
    bl.means.sims <- matrix(ncol=length(deme), nrow=length(simNum))
    bl.vars.sims <- matrix(ncol=length(deme), nrow=length(simNum))
    #junction density is stored with a row for each simulation and junction deme combination
    # and each row needs the deme number and samplesize appended to the beginning
    junden.over.sims <- matrix(nrow=length(simNum)*length(junction_demes), ncol=chrLength[chr+1]/winSize + 2 )
    hapblocks.m.over.sims <- matrix(nrow = length(simNum)*length(junction_demes), ncol= sampSize * 2 + 2 )
    hapblocks.v.over.sims <- matrix(nrow = length(simNum)*length(junction_demes), ncol= sampSize * 2 + 2 )
    simCount <- 0
    for(s in simNum)
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
      
      #
      
      # calculate the mean and variance of block length for the entire metapopulation
      meta.mean_bl <- mean(blockLengths$Length)
      meta.var_bl <- var 
      
      # calculate the mean and variance in block length for each deme individually
      deme.mean_bl <- c()
      deme.var_bl <- c()
      for(d in deme)
      {
        thisDeme <- subset(blockLengths, blockLengths$Deme==d) # get blocks only from the deme we are currently analyzing
        thisDeme <- subset(thisDeme, thisDeme$Type!=3) # exclude any Y-chromosomes
        thisDeme <- na.omit(thisDeme)
        bl.means.sims[s+1,d+1] <- mean(thisDeme$Length)
        bl.vars.sims[s+1,d+1] <- var(thisDeme$Length)
      }
      
      # calculate the junction density in windows across the genome
      # first make sure the window size is not bigger than the actual chromosome
      if(winSize > chrLength[chr+1]) # if it is
      {
        windows <- c(chrLength[chr+1]) # we only need one window
      }
      else # if not
      {
        windows <- seq(winSize, chrLength[chr+1], winSize) # then make list of all the windows by their end position
      }
      
      demeCount <- 0
      for( jd in junction_demes )
      {
        index <- s * length(junction_demes) + jd - junction_demes[1] + 1
        meta.junction_density <- rep(0,length(windows))
        
        thisDeme.bl <- subset(blockLengths, blockLengths$Deme == jd ) # use only the current deme
        thisDeme.bl <- subset(thisDeme.bl, thisDeme.bl$Type != 3 ) # and exclude any Y-chromosomes
        thisDeme.bl <- na.omit(thisDeme.bl)
        thisDeme.hap_means.bl <- c(s,jd,rep(0,sampSize*2))
        thisDeme.hap_vars.bl <- c(s,jd,rep(0,sampSize*2))
        hapBlockIndex <- simCount * length(junction_demes) + demeCount
        
        blocks.by.hap <- c()
        lastHap <- thisDeme.bl[1,3] # assign the first block seen as the previous block
        sampleSize <- 1 # start counting the number of chromosomes that are analyzed
        
        # Go through all the blocks in this deme and count the junctions in each window
        for(h in 1:nrow(thisDeme.bl) )
        {
          # check to see if the current block/junction is at the start of the chromosome
          if(thisDeme.bl[h,7] != 0 ) # if it is not
          {
            win <- floor(thisDeme.bl[h,7]/winSize +1 ) # calculate the index of the window
            meta.junction_density[win] <- meta.junction_density[win] + 1 # and add 1 to the count for that window
          }
          
          # Now check to see if this block is on the same haplotype as the last one
          thisHap <- thisDeme.bl[h,3]
          if(thisHap == lastHap ) # if they match
          {
            blocks.by.hap <- c(blocks.by.hap, thisDeme.bl[h,11]) # add its length to the list for that haplotype
            lastHap <- thisHap # advance the lastHap to the current without doing anything else
          }
          else
          {
            thisDeme.hap_means.bl[2+sampleSize] <- mean(blocks.by.hap) 
            thisDeme.hap_vars.bl[2+sampleSize] <- var(blocks.by.hap)
            lastHap <- thisHap
            sampleSize <- sampleSize + 1
          }
          
        }
        
        # convert each junction density to the average over haplotypes for the deme
        for( i in 1:length(meta.junction_density) )
        {
          meta.junction_density[i] <- meta.junction_density[i] / sampleSize
        }
        
        junden.over.sims[index,] <- c(jd, sampleSize, meta.junction_density )
        hapblocks.m.over.sims[hapBlockIndex,] <- thisDeme.hap_means.bl
        hapblocks.v.over.sims[hapBlockIndex,] <- thisDeme.hap_vars.bl
        thisDeme.hap_means.bl <- c()
        thisDeme.hap_vars.bl <- c()
        demeCount <- demeCount + 1
      }
      simCount <- simCount + 1
    }
    
    # Save all the matrices with the data to files
    save(junden.over.sims, file=paste(expName, "_", jobID, "/", expName, "_", chr, "_", gen, "_jd.R", sep=""))
    save(bl.means.sims, file=paste(expName, "_", jobID, "/", expName, "_", chr, "_", gen, "_blm.R", sep=""))
    save(bl.vars.sims, file=paste(expName, "_", jobID, "/", expName, "_", chr, "_", gen, "_vlm.R", sep=""))
    save(hapblocks.m.over.sims, file=paste(expName, "_", jobID, "/", expName, "_", chr, "_", gen, "_hbm.R", sep=""))
    save(hapblocks.v.over.sims, file=paste(expName, "_", jobID, "/", expName, "_", chr, "_", gen, "_hbv.R", sep=""))
    
  }
}
