library(foreach)
library(iterators)
library(doParallel)

#Register the parallel backend
registerDoParallel(16)

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/h23v14/')

dir.create(paste(getwd(),'/DEM/',sep=''))
dir.create(paste(getwd(),'/IMG/',sep=''))
dir.create(paste(getwd(),'/MAPS/',sep=''))
dir.create(paste(getwd(),'/NLCD/',sep=''))
dir.create(paste(getwd(),'/PHENO/',sep=''))
dir.create(paste(getwd(),'/SHP/',sep=''))

in_dirs_SR <- list.files(path=getwd(),pattern=glob2rx("*SR*"),
  full.names=T,include.dirs=T)
in_dirs_TA <- list.files(path=getwd(),pattern=glob2rx("*TA*"),
  full.names=T,include.dirs=T)

pred.SPR <- foreach(i = 1:length(in_dirs_TA), .combine = rbind) %dopar% {
  print(i)

  dir_name <- substr(in_dirs_TA[i],57,88)
  dir.create(paste(getwd(),'/IMG/',dir_name,sep=''))

  untar(in_dirs_SR[i],exdir = paste(getwd(),'/IMG/',dir_name,sep=''))
  untar(in_dirs_TA[i],exdir = paste(getwd(),'/IMG/',dir_name,sep=''))

  unlink(in_dirs_SR[i])
  unlink(in_dirs_TA[i])

  in_dirs_TAB <- list.files(path=paste(getwd(),'/IMG/',dir_name,sep=''),
    pattern=glob2rx("*TAB*"),full.names=T,include.dirs=T)
  unlink(in_dirs_TAB)

  pSPR <- i
}


