library(rgdal)
library(raster)

img <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/h22v15/LC08_CU_022015_20130415_20170803/evi2.tif')

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/h22v15/PHENO')
tmp_files <- list.files(pattern = "evi2_phenology*", recursive = TRUE,full.names = TRUE)
chunk <- unlist(lapply(tmp_files, 
  function(x) na.omit(as.numeric(unlist(strsplit(unlist(x), "[^0-9]+"))))[2]))

phen <- matrix(NA,ncell(img),4)
spr_anom <- matrix(NA,ncell(img),1)

for (j in 1:length(tmp_files)){
  print(j)
  load(tmp_files[j])
  
  phen[(25*5000*(chunk[j]-1)+1):((25*5000)*chunk[j]),] <- all_pheno[,1:4]
  
  spr_anom[(25*5000*(chunk[j]-1)+1):((25*5000)*chunk[j]),] <- all_pheno[,35]-all_pheno[,3]
}

phen[phen<=0] <- NA
phen[which(phen[,2]<0.95),3:4] <- NA
spr_anom[which(phen[,2]<0.95),] <- NA

phen[,2] <- round(phen[,2]*100)

names <- c('nobs','rsquare','spr','aut')
for (k in 1:4){
  print(k)
  s <- setValues(img,phen[,k])
  name <- paste(names[k],'.tif',sep="")
  writeRaster(s,filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/h22v15/MAPS/',
    name,sep=""),format='GTiff',overwrite=TRUE)
}

names <- c('spr2012anom')
for (k in 1:1){
  print(k)
  s <- setValues(img,spr_anom[,k])
  name <- paste(names[k],'.tif',sep="")
  writeRaster(s,filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/h22v15/MAPS/',
    name,sep=""),format='GTiff',overwrite=TRUE)
}