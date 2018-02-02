# Get the distribution of junction densities for a locus on a chromosome
# input = c(chromosome, position, generation, deme )
junden_dist_locus <- function(y)
{
  x <- y
  for( i in 1:1000)
  {
    x <- c(x, sim_jd[[i]][[y[3]]][[y[1]]][[y[4],y[2]+2]])
  }
  return(x)
}


#calculate the mean number of junctions
calc_jd_sim_mean <- function(s, g, c, d, w) # sim(s), deme(S), chromosome(s), window(s)
{
  w <- w + 2 # the first two positions are chr, and postion
  jd_mean <- 0
  count <- 1
  
  for( i in s )
  {
    for( j in g )
    {
      for ( k in c )
      {
        for ( l in w )
        {
          jd_mean <- jd_mean + sim_jd[[i]][[j]][[k]][d,l]
          #print(paste(i,j,k,l,jd_mean,sep=" "))
          count <- count + 1
        }
      }
    }
  }
  jd_mean <- jd_mean / count
  return(jd_mean)
}

#calculate the variance AMONG JD WINDOWS in the number of junctions
calc_jd_sim_var <- function(s, g, c, d, w) # sim(s), deme(S), chromosome(s), window(s)
{
  w <- w + 2 # the first two positions are chr, and postion
  reserve <- length(s) * length(g) * length(c) * length(w)
  jd_var <- rep(0,reserve)
  count <- 1
  #print(reserve)
  
  for( i in s )
  {
    for( j in g )
    {
      for ( k in c )
      {
        for ( l in w )
        {
          jd_var[count] <- sim_jd[[i]][[j]][[k]][d,l]
          #print(paste(i,j,k,l,jd_var[count],sep=" "))
          count <- count + 1
        }
      }
    }
  }
  #print(jd_var)
  return(var(jd_var))
}

#calculate the Confidence Interval AMONG JD WINDOWS in the number of junctions
calc_jd_sim_CI <- function(s, g, c, d, w, ci) # sim(s), deme(S), chromosome(s), window(s)
{
  w <- w + 2 # the first two positions are chr, and postion
  reserve <- length(s) * length(g) * length(c) * length(w)
  jd_CI <- rep(0,reserve)
  count <- 1
  
  for( i in s )
  {
    for( j in g )
    {
      for ( k in c )
      {
        for ( l in w )
        {
          jd_CI[count] <- sim_jd[[i]][[j]][[k]][d,l]
          #print(paste(i,j,k,l,jd_CI[count],sep=" "))
          count <- count + 1
        }
      }
    }
  }
  
  jd_CI <- sort(jd_CI)
  
  top <- round( length(jd_CI) - ( 0.5 * (1-ci ) * length(jd_CI) )  )
  bottom <- round( 0.5 * (1 - ci) * length(jd_CI) )
  
  #print(bottom,top,sep=" ")
  
  if(bottom < 1 ) 
  { 
    bottom <- 1
  }
  
  return( c( jd_CI[bottom], jd_CI[top] ) )
}

#FUNCTION: calculate the mean block length and 95% CI for a given chromosome, generation, and deme
calc_bl_sim_mCI <- function(s, c, g, d, ci)
{
  m <- c()
  for(s in sims)
  {
    n <- sim_mbl[[s]][[g]][[c]][[d]]
    if(c==1)
    {
      
    }
    m <- c(m, n)
  }
  
  m <- sort(m)
  top <- round( length(m) - (1-ci)/2 * length(m))
  
  bottom <- round( (1-ci)/2 * length(m) )
  
  if(bottom < 1)
  {
    bottom <- 1
  }
  
  return( c(mean(m),m[bottom],m[top] )  )
}

#FUNCTION: calculate the mean block length and 95% CI for a given chromosome, generation, and deme
calc_bl_sim_vCI <- function(s, c, g, d, ci)
{
  m <- c()
  for(s in sims)
  {
    n <- sim_vbl[[s]][[g]][[c]][[d]]
    m <- c(m, n)
  }
  
  m <- sort(m)
  top <- round( length(m) - (1-ci)/2 * length(m))
  
  bottom <- round( (1-ci)/2 * length(m) )
  
  if(bottom < 1)
  {
    bottom <- 1
  }
  
  return( c(mean(m),m[bottom],m[top] )  )
}

#make a data table, to be plotted, of mean junction density with confidence intervals
num_sims <- 1000
sims <- passed_sims + 1
windows <- seq(1,10,1)
chromosome <- c(2)
deme <- c(5)
generation <- c(11)
mean_jd_table <- matrix(nrow=1,ncol=7)
colnames(mean_jd_table) <- c("Chr","Pos","Gen","Deme","Mean_jd","Lower","Upper")
for(g in generation)
{
  for(d in deme)
  {
    for(chr in chromosome)
    {
      for(win in windows)
      {
          new_row <- matrix(nrow=1,ncol=7)
          new_row[,1] <- chr
          new_row[,2] <- win / length(windows)
          new_row[,3] <- g
          new_row[,4] <- d
          new_row[,5] <- calc_jd_sim_mean(sims,generation,chr,d,win)
          new_row[,c(6,7)] <- calc_jd_sim_CI(sims,generation,chr,d,win,0.95)
          mean_jd_table <- rbind(mean_jd_table, new_row)
      }
    }
  }
}

mean_jd_table <- mean_jd_table[-1,]
#colnames(mean_jd_table) <- c("Chr","Pos","Gen","Deme","MeanJD","Lower","Upper")
mean_jd_table <- data.frame(mean_jd_table)
save(mean_jd_table, file=paste(expName,"_jdtab.R",sep="") )

#make a data table, to be plotted, of mean and variance of block length with confidence intervals
generations <- c(2:13)
mean_bl <- matrix(nrow=1,ncol=9)
colnames(mean_bl) <- c("Chr","Gen","Deme","Mean_bl","L_MBL","U_MBL","Var_bl","L_VBL","U_VBL")
for(d in deme)
{
  for(c in chromosome)
  {
    for(g in generations)
    {
      new_row <- matrix(nrow=1, ncol=9)
      new_row[,1] <- c
      new_row[,2] <- g
      new_row[,3] <- d
      new_row[,c(4,5,6)] <- calc_bl_sim_mCI(sims, c, g, d, 0.95)
      new_row[,c(7,8,9)] <- calc_bl_sim_vCI(sims, c, g, d, 0.95)
      mean_bl <- rbind(mean_bl, new_row)
    }
  }
}

mean_bl <- mean_bl[-1,]
#colnames(mean_jd_table) <- c("Chr","Pos","Gen","Deme","Mean_jd","Lower","Upper")
mean_bl <- data.frame(mean_bl)
save(mean_bl, file=paste(expName,"_bltab.R",sep="") )
