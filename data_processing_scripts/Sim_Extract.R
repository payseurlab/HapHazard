
 
expNames <-  c( "Jon_S" , "Benjen_S")
  
JDWS_list <- list( c( "Jon_S" = c(), "Benjen_S" = c() ) )

# Get the distribution of junction densities for a locus on a chromosome
junden_dist_locus <- function(y)
{
  x <- y
  for( i in 1:1000)
  {
    x <- c(x, sim_jd[[i]][[y[3]]][[y[1]]][[y[4],y[2]+2]])
  }
  return(x)
}

# calculate a confidence interval on a vector x whose range is given by ci
calc_CI <- function(x, ci)
{
  x <- sort(x)
  y <- c()
  
  top <- round( length(x) - (1-ci)/2 * length(x))
  
  bottom <- round( (1-ci)/2 * length(x) )
  
  if(bottom < 1)
  {
    bottom <- 1
  }
  
  y <- c(x[bottom],x[top])
  
  return(y)
}
  
# jd_win input = c(chromosome, position, generation, deme )
jd_win1 <- c(1,50,10,5)
jd_win2 <- c(2,50,10,5)
jd_win3 <- c(1,50,10,6)
jd_win4 <- c(2,50,10,6)

for(expName in expNames)
{
num_of_sims <- 1000 - 1
sims <- seq(0,num_of_sims,1)
tar_command <- paste("tar -xzf ",expName,".tar.gz",sep="")
system(tar_command)

#clines, age, freq spec, mbl, vbl, junden, hap-mbl, hap-vbl, hyb index
data_toggle <- c(FALSE,FALSE,FALSE,TRUE,TRUE, TRUE, TRUE, TRUE, TRUE)

sim_cline <- list()
sim_age <- list()
sim_fs <- list()
sim_mbl <- list()
sim_vbl <- list()
sim_jd <- list()
sim_hbm <- list()
sim_hbv <- list()
sim_hyb <- list()

passed_sims <- c()
failed_sims <- c()
failed_tars <- c()
missing_files <- c()


for(s in sims)
{
  dir <- paste("~/Desktop/Theo_Desktop/",expName,"/",expName,"_",s,"/",sep="")
  print(dir)
  setwd(dir)
  
  tar_file <- paste(expName,"_",s,".tar.gz",sep="")
  
  if( file.exists(tar_file) )
  {
    passed_sims <- c(passed_sims, s)
    tar_command <- paste("tar -xzf ",expName,"_",s,".tar.gz",sep="")
    system(tar_command)
  
    dir <- paste(dir,expName,"_",s,"/",sep="")
    setwd(dir)
    source("Rvars.R")
    
    generation_cline <- list()
    generation_age <- list()  
    generation_fs <- list()
    generation_mbl <- list()
    generation_vbl <- list()
    generation_jd <- list()
    generation_hbm <- list()
    generation_hbv <- list()
    generation_hyb <- list()
    g_index <- 1
    
    for(g in generations )
    {
      cline_file <- paste(g,"bfc.R",sep=".")
      if(file.exists(cline_file) && data_toggle[1] )
      {
        load(cline_file)
        generation_cline[[g_index]] <- cline.bestfit
      }
      else
      {
        if(!file.exists(cline_file))
        {
          missing_files <- c(missing_files, cline_file)
        }
        generation_cline <- NA
      }
      
      chromosome_age <- list()
      chromosome_fs <- list()
      chromosome_mbl <- list()
      chromosome_vbl <- list()
      chromosome_jd <- list()
      chromosome_hbm <- list()
      chromosome_hbv <- list()
      chromosome_hyb <- list()
      
      for( c in chromosomes )
      {
        age_file <- paste(expName,"_0_",c,"_",g,"_age.R",sep="")
        
        if(file.exists(age_file) && data_toggle[2] )
        {
          load(age_file)
          chromosome_age[[1+c]] <- Age.by.deme
        }
        else
        {
          if(!file.exists(age_file))
          {
            missing_files <- c(missing_files, age_file)
          }
          chromosome_age[[1+c]] <- NA
        }
        
        fs_file <- paste(expName,"_0_",c,"_",g,"_frs.R",sep="")
        if(file.exists(fs_file) & data_toggle[3] )
        {
          load(fs_file)
          chromosome_fs[[1+c]] <- Freq.Spec
        }
        else
        {
          if(!file.exists(fs_file))
          {
            missing_files <- c(missing_files, fs_file)
          }
          chromosome_fs[[1+c]] <- NA
        }
        
        mbl_file <- paste(expName,"_",c,"_",g,"_blm.R",sep="")
        if(file.exists(mbl_file) && data_toggle[4] )
        {
          load(mbl_file)
          chromosome_mbl[[1+c]] <- bl.means.sims
        } 
        else
        {
          if(!file.exists(mbl_file))
          {
            missing_files <- c(missing_files, mbl_file)
          }
          chromosome_mbl[[1+c]] <- NA
        }
        
        vbl_file <- paste(expName,"_",c,"_",g,"_vlm.R",sep="")
        if(file.exists(vbl_file) && data_toggle[5] )
        {
          load(vbl_file)
          chromosome_vbl[[1+c]] <- bl.vars.sims
        }
        else
        {
          if(!file.exists(vbl_file))
          {
            missing_files <- c(missing_files, vbl_file)
          }
          chromosome_vbl[[1+c]] <- NA
        }
        
        jd_file <- paste(expName,"_",c,"_",g,"_jd.R",sep="")
        if(file.exists(jd_file) && data_toggle[6])
        {
          load(jd_file)
          chromosome_jd[[1+c]] <- junden.over.sims          
        }
        else
        {
          if(!file.exists(jd_file))
          {
            missing_files <- c(missing_files, jd_file)
          }
          chromosome_jd[[1+c]]
        }
          
        hbm_file <- paste(expName,"_",c,"_",g,"_hbm.R",sep="")
        if(file.exists(hbm_file) && data_toggle[7])
        {
          load(hbm_file)
          chromosome_hbm[[1+c]] <- hapblocks.m.over.sims
        }
        else
        {
          if(!file.exists(hbm_file))
          {
            missing_files <- c(missing_files, hbm_file)
          }
          chromosome_hbm[[1+c]]
        }
        
        hbv_file <- paste(expName,"_",c,"_",g,"_hbv.R",sep="")
        if(file.exists(hbv_file)&& data_toggle[8])
        {
          load(hbv_file)
          chromosome_hbv[[1+c]] <- hapblocks.v.over.sims
        }
        else
        {
          if(!file.exists(hbv_file))
          {
            missing_files <- c(missing_files, hbv_file)
          }
          chromosome_hbv[[1+c]] <- NA
        }
        
        hyb_file <- paste("BL_CHR",c,"/0.",g,".hyb",sep="")
        if(file.exists(hyb_file)&& data_toggle[9])
        {
          chromosome_hyb[[1+c]] <-  read.csv(hyb_file)
        }
        else
        {
          if(!file.exists(hyb_file))
          {
            missing_files <- c(missing_files, hyb_file)
          }
          chromosome_hyb[[1+c]] <- NA
        }
      }
      generation_age[[g_index]] <- chromosome_age
      generation_fs[[g_index]] <- chromosome_fs
      generation_mbl[[g_index]] <- chromosome_mbl
      generation_vbl[[g_index]] <- chromosome_vbl
      generation_jd[[g_index]] <- chromosome_jd
      generation_hbm[[g_index]] <- chromosome_hbm
      generation_hbv[[g_index]] <- chromosome_hbv
      generation_hyb[[1+c]] <- chromosome_hyb
      g_index <- g_index + 1
    }
    
    sim_cline[[s+1]] <- generation_cline
    sim_age[[s+1]] <- generation_age
    sim_fs[[s+1]] <- generation_fs
    sim_mbl[[s+1]] <- generation_mbl
    sim_vbl[[s+1]] <- generation_vbl
    sim_jd[[s+1]] <- generation_jd
    sim_hbm[[s+1]] <- generation_hbm
    sim_hbv[[s+1]] <- generation_hbv
    sim_hyb[[s+1]] <- generation_hyb
  }
  else
  {
    print(paste("Could not open sim ",tar_file,",tar file did not exist.",sep=""))
    failed_tars <- c(failed_tars, tar_file)
    failed_sims <- c(failed_sims, s)
  }
}
setwd("~/Desktop/Theo_Desktop")

#source("Sim_Summary.R")

jd_win_sum <- junden_dist_locus(jd_win1)
jd_win_sum <- rbind(jd_win_sum, junden_dist_locus(jd_win2))
jd_win_sum <- rbind(jd_win_sum, junden_dist_locus(jd_win3))
jd_win_sum <- rbind(jd_win_sum, junden_dist_locus(jd_win4))

JDWS_list[[expName]] <- jd_win_sum

save(jd_win_sum, file=paste(expName,"_jdws.R",sep=""))


}
#clean up the workspace
#rm(bl.means.sims, bl.vars.sims,cline.bestfit,hapblocks.m.over.sims,hapblocks.v.over.sims,junden.over.sims)
#rm(age_file,c,chrLength,chromosome_age,chromosome_fs,chromosome_hbm,chromosome_hbv,chromosome_hyb,chromosome_jd,chromosome_mbl,chromosome_vbl)
#rm(generation_cline,generation_age,generation_fs,generation_hbm,generation_hbv,generation_hyb,generation_jd,generation_mbl,generation_vbl)
#rm(chromosomes,cline_file,data_toggle,demes,dir,g,g_index,generations)
#rm(fs_file,hbm_file,hbv_file,hyb_file,jd_file,jobID,junction_demes,mbl_file,s,num_of_sims,sims,tar_command)
#rm(tar_file,vbl_file,winSize)