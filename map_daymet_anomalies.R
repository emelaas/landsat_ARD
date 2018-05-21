library(rgdal)
library(raster)


my.colors3 = colorRampPalette(c("firebrick4","orangered","orange","white",
  "cyan","deepskyblue","deepskyblue4"))
my.colors4 = colorRampPalette(c("deepskyblue4","deepskyblue","cyan","white",
  "orange","orangered","firebrick4"))

mod <- c('Chill','Photo','SW_Jan1','SW_Mar13')

LTM <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sprLTM.tif')

# Generate LTM prediction map for each model
setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
for (m in 1:4){
  print(m)

  in_dirs_tile <- list.files(path=getwd(),
    pattern=glob2rx(paste("spr_mod",m,"*tif",sep="")),
    full.names=T,include.dirs=T,recursive=TRUE)
  all_mod <- stack(in_dirs_tile)
  all_mod_vals <- getValues(all_mod)
  all_mod_mean <- rowMeans(all_mod_vals,na.rm=TRUE)
  mod_mean <- setValues(LTM,all_mod_mean)
  writeRaster(mod_mean,
    filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr_mod',m,'_LTM',sep=""),
    format='GTiff',overwrite=TRUE)
}


for (year in 1984:2017){

  print(year)

  #x11(h=8,w=11)
  pdf(h=8,w=11,
    paste('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/spr_mod2_1km_anomaly_',
    year,'.pdf',sep=""))
  par(mfrow=c(2,3),mar=c(3,3,3,3))

  LTM <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sprLTM.tif')
  y_obs <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr',year,'.tif',sep=""))
  diff <- y_obs-LTM
  diff[diff < -21] <- -21
  diff[diff > 21] <- 21
  plot(diff,main=paste('Obs. Anomaly ',year,sep=''),colNA='black',col=my.colors3(150),
    zlim=c(-21,21),axes=FALSE,bty='n',legend.width=1.2)

  Zscore <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/agdd_Zscore_',year,'.tif',sep=""))
  Zscore[is.na(LTM)==1] <- NA
  plot(Zscore,main='AGDD Z Score',colNA='black',col=my.colors4(150),
    zlim=c(-3,3),axes=FALSE,bty='n',legend.width=1.2)

  for (m in 1:4){

    print(m)

    LTM_mod <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr_mod',m,'_LTM.tif',sep=""))

    y_pred <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr_mod',m,
      '_',year,'.tif',sep=""))
    diff <- y_pred-LTM_mod
    diff[diff < -21] <- -21
    diff[diff > 21] <- 21

    diff[is.na(LTM)==1] <- NA

    plot(diff,main=mod[m],colNA='black',col=my.colors3(150),
      zlim=c(-21,21),axes=FALSE,bty='n',legend.width=1.2)
  }

  dev.off()
}
