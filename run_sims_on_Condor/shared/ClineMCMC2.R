# The name of the experiment to be analyzed is given on the command line
#args <- commandArgs(TRUE)

#for(i in 1:length(args))
#{
#	eval(parse(text=args[[i]]))
#}
#demes <- c(0:9)
#expName <- "Benjen_S"
#generations <- c(499,999,1499,1999)

print("Running ClineMCMC2.R")

# Logistic cline -- calculate the logistic function over x
logit <- function(x, w, c)
{
  cline <- 1 / (1 + exp(-w*(x - c)))
  return(cline)
}

# select a random value of x from a normal distribution with mean m, std dev s, and require the variable be bound by
# upper and lower limits
p_norm_bound <- function(m, s, lower, upper)
{
  x <- rnorm(1,mean=m,sd=s)
  
  while( ( x < lower ) || ( x > upper ) )
  {
    #print(paste(x,lower,upper,sep=" "))
    x <- rnorm(1,mean=m,sd=s)
  }
  
  return(x)
}

for(gen in generations )
{
  
  filename <- paste(expName,"_",jobID,"/", expName, ".0.",gen,".clines",sep="")
  filename2 <- paste(expName,"_",jobID,"/",expName, ".0.",gen,".cline0",sep="")
  filename3 <- paste(expName,"_",jobID,"/",expName, ".0.",gen,".clines0",sep="")

  if(file.exists(filename))
  {
    cline.data <- read.csv(filename) 
  }
  else if (file.exists(filename2))
  {
    cline.data <- read.csv(filename2)
  }
  else if (file.exists(filename3))
  {
    cline.data <- read.csv(filename3)
  }
  else
  {
    print(paste("Could not open cline file, ",filename,"!",sep=""))
    break;
  }

  m <- nrow(cline.data)
  n <- 5
  cline.bestfit <- matrix(rep(100, m*n), nrow=m, ncol=n)
  
  center <- 4.5
  width <- 1
  
  for(j in 1:m)
  {
    x <- cline.data[j,]
    chr <- as.numeric(x[1])
    pos <- as.numeric(x[2])
    x <- x[-c(1:2)]
    y <- ( 2 * demes + 2) 
    x <- x[,y] / 20
    
    burn <- c(0:1000)
    
	if( length(demes) < 2 )
	{
		cline.bestfit[j,] <- c( chr, pos, 0, 0, 0 )
	}
	else
	{

		for( i in burn )
		{
		  center <- p_norm_bound(center, 1, 0, 9)
		  width <- p_norm_bound(width, 1, 0, 10)
		  cline <- logit(demes,width,center)
		  rmse <- sqrt( mean( (x - cline)^2 ) )
		  
		  if( rmse < cline.bestfit[j,5] )
		  {
		    cline.bestfit[j,] <- c(chr, pos, width, center, rmse)
		  }
		}
	}
  }
  colnames(cline.bestfit) <- c("Chr","Pos","Width","Center","RMSE")
  save(cline.bestfit, file=paste(expName,"_",jobID,"/",gen,".bfc.R",sep=""))

}
