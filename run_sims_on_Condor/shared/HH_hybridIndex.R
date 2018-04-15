
# The name of the experiment to be analyzed is given on the command line
#args <- commandArgs(TRUE)

#for(i in 1:length(args))
#{
	#eval(parse(text=args[[i]]))
#}

#first open the block length/junction file

print("Running HH_hybridIndex.R")

#chromosomes <- seq(0,args[2],1)
#generations <- c(249,499,749,999,1499,1999)
#deme <- c(0:9)
#winSize <- 0.01 # in Morgans
center_demes <- demes

for(chr in chromosomes)
{
  for(gen in generations)
  {
    list.index <- 4
    markerFileName <- paste(expName, "_",jobID,"/BL_CHR",chr,"/0.",gen,".mkr", sep="")

    if(file.exists(markerFileName))
    {
      print(markerFileName)
      markers <- read.csv(markerFileName)
      markers <- markers[-nrow(markers),]
    }
    else
    {    
      print(paste("Could not open file: ", markerFileName, sep=""))
    }
    
    r <- 1
    
    hets_HI <- rep.int(0,nrow(markers))
    while ( r <= nrow(markers) )
    {
      q <- 6
      het <- c(0,0)
      hyb <- c(0,0)
      while( q <= ncol(markers) )
      {
        if( markers[r,5] != 3 && markers[r+1,5] != 3 )
        {
          if( markers[r,q] != markers[r+1,q] )
          {
            het[1] <- het[1] + 1
          }
          het[2] <- het[2] + 1
          
          hyb[1] <- hyb[1] + markers[r,q] + markers[r+1,q]
          hyb[2] <- hyb[2] + 2
        }
        else
        {
          if( markers[r,q] == 3 )
          {
            hyb[1] <- hyb[1] + markers[r+1,q]
          }
          else
          {
            hyb[1] <- hyb[1] + markers[r,q]
          }
          hyb[2] <- hyb[2] + 1
        }
        q <- q + 1
      }
      
      hets_HI[r] <- het[1]/het[2]
      hets_HI[r+1] <- hyb[1]/hyb[2] 
      r <- r + 2
      
    }

    hetz.HI <- markers[seq(1,nrow(markers),2),c(1:4)]
    heterozygosity <- hets_HI[seq(1,length(hets_HI),2)]
    hybridIndex <- hets_HI[seq(2,length(hets_HI),2)]
    hetz.HI <- cbind(hetz.HI, heterozygosity, hybridIndex)
    
    write.csv(hetz.HI, file=paste(expName, "_",jobID,"/BL_CHR",chr,"/0.",gen,".hyb", sep=""),quote=FALSE,row.names=FALSE)
    
  }
  
}



