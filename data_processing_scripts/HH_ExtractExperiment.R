whole_tar_filename <- paste(experiment_name,'.tar.gz',sep="")
if ( file.exists(whole_tar_filename)) {
  whole_tar_command <- paste("tar -xzf ",whole_tar_filename,sep="")
  system(whole_tar_command)
} else {
  print(paste("No",whole_tar_filename,"found"))
}



sims_not_found <- c()

for(sim in 1:number_of_sims)
{
  exp_dir <- paste(workdir,experiment_name,"/",experiment_name,"_",sim,"/",sep="")
  if( dir.exists(exp_dir))
  {
    setwd(exp_dir)
    tar_filename <- paste(experiment_name,"_",sim,'.tar.gz',sep="")
    if ( file.exists(tar_filename))
    {
      tar_command <- paste("tar -xzf ", tar_filename,sep="")
      system(tar_command)
    }
    else
    {
      not_found <- paste(exp_dir,experiment_name,'_',sim,'.tar.gz',sep="")
      sims_not_found <- c(sims_not_found, not_found)
    }
  }
  else
  {
    not_found <- paste(exp_dir,experiment_name,'_',sim,'.tar.gz',sep="")
    sims_not_found <- c(sims_not_found, not_found)
  }
}

setwd(workdir)

