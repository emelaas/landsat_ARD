library(rgdal)
library(raster)

can_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/Canada/','Canada')
can_shp_proj <- spTransform(can_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
usa_shp_proj <- spTransform(usa_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_crop <- crop(usa_shp_proj,extent(194214,2426712,-1660235,946090))
can_crop <- crop(can_shp_proj,extent(194214,2426712,-1660235,946090))

my.colors3 = colorRampPalette(c("firebrick4","orangered","orange","white",
  "cyan","deepskyblue","deepskyblue4"))
my.colors4 = colorRampPalette(c("deepskyblue4","deepskyblue","cyan","white",
  "orange","orangered","firebrick4"))

mod <- c('Chill','Photo','SW Mar 17','SW Jan 1 (Mean AGDD)')

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
    paste('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/model_bias_by_yr/spr_mod2_1km_bias_',
    year,'.pdf',sep=""))
  par(mfrow=c(2,3),mar=c(3,3,3,3))

  LTM <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sprLTM.tif')
  y_obs <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr',year,'.tif',sep=""))
  diff <- y_obs-LTM
  diff[diff < -14] <- -14
  diff[diff > 14] <- 14
  
  plot(usa_crop,col='white',main=paste('Obs. Anomaly ',year,sep=''))
  rect(194214,-1660235,2426712,946090,col='lightblue')
  plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(diff,col=my.colors3(150),zlim=c(-14,14),axes=FALSE,bty='n',legend.width=1.2,add=TRUE)

  Zscore <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/agdd_Zscore_',year,'.tif',sep=""))
  Zscore[is.na(LTM)==1] <- NA
  plot(usa_crop,col='white',main='AGDD Z Score')
  rect(194214,-1660235,2426712,946090,col='lightblue')
  plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(Zscore,col=my.colors4(150),zlim=c(-3,3),axes=FALSE,bty='n',legend.width=1.2,add=TRUE)

  for (m in 1:4){

    print(m)

    LTM_mod <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr_mod',m,'_LTM.tif',sep=""))

    y_pred <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr_mod',m,
      '_',year,'.tif',sep=""))
    diff <- y_pred-y_obs
    diff[diff < -14] <- -14
    diff[diff > 14] <- 14

    diff[is.na(LTM)==1] <- NA

    plot(usa_crop,col='white',main=mod[m])
    rect(194214,-1660235,2426712,946090,col='lightblue')
    plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
    plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
    plot(diff,col=my.colors3(150),zlim=c(-14,14),axes=FALSE,bty='n',legend.width=1.2,add=TRUE)
  }

  dev.off()
}


#Plot density of bias between each model and annual observed phenology from Landsat
# - only include pixels with AGDD Z scores > 2
below5 <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/below5_sum.tif')
below5_vals <- getValues(below5) 

# x11(h=6.5,w=13)
pdf(h=6.5,w=13,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/model_bias_by_yr/density_extreme_years.pdf')
par(mfrow=c(2,5),mar=c(4,3,2,3))

for (year in c(1987,1998,2010,2012,2017)){
  Zscore <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/agdd_Zscore_',year,'.tif',sep=""))
  Zscore_vals <- getValues(Zscore)
  w <- which(Zscore_vals>2)
  
  print(year)
  
  y_obs <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr',year,'.tif',sep=""))
  
  setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
  in_dirs_ypred <- list.files(path=getwd(),
    pattern=glob2rx(paste("spr_mod*",year,"*tif",sep='')),full.names=T,include.dirs=T,recursive=TRUE)
  
  y_pred <- stack(in_dirs_ypred)
  s <- stack(replicate(4, y_obs))
  diff <- s-y_pred
  diff_vals <- getValues(diff)
  colnames(diff_vals) <- c('Chill','Photo','SW Mar 17','SW Jan 1 (Mean AGDD)')
  
  diff_vals1 <- diff_vals[which(Zscore_vals>2),]
  
  if (length(diff_vals1)>4){
    dens <- apply(diff_vals1, 2, density, bw=2, na.rm=TRUE)
    plot(NA, xlim=c(-21,21), ylim=range(sapply(dens, "[", "y")),
      main=year,
      xlab='Model Bias (days)')
    mapply(lines, dens, col=1:length(dens),lwd=4)
    abline(v=0,lty=2)
  }
  
  if (year==2017) legend("topleft", legend=names(dens), fill=1:length(dens),bty='n')

}

for (year in c(1987,1998,2010,2012,2017)){
  Zscore <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/agdd_Zscore_',year,'.tif',sep=""))
  plot(usa_crop,col='white',main='')
  rect(194214,-1660235,2426712,946090,col='lightblue')
  plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(Zscore,col=my.colors4(150),zlim=c(-3,3),axes=FALSE,bty='n',legend.width=1.2,add=TRUE)
}

dev.off()

for (year in 1984:2017){
  print(year)
  
  y_obs <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr',year,'.tif',sep=""))
  
  setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
  in_dirs_ypred <- list.files(path=getwd(),
    pattern=glob2rx(paste("spr_mod*",year,"*tif",sep='')),full.names=T,include.dirs=T,recursive=TRUE)
  
  y_pred <- stack(in_dirs_ypred)
  s <- stack(replicate(4, y_obs))
  diff <- y_pred-s
  diff_vals <- getValues(diff)
  colnames(diff_vals) <- c('Chill','Photo','SW Mar 17','SW Jan 1 (Mean AGDD)')
  
  Zscore <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/agdd_Zscore_',year,'.tif',sep=""))
  Zscore_vals <- getValues(Zscore)
  
  diff_vals1 <- diff_vals[which(Zscore_vals<1),]
  diff_vals2 <- diff_vals[which(Zscore_vals>=1 & Zscore_vals<2),]
  diff_vals3 <- diff_vals[which(Zscore_vals>2),]
  
  x11(h=4,w=11)
#   pdf(h=4,w=11,
#     paste('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/model_bias_by_yr/spr_mod2_1km_bias_density',
#       year,'.pdf',sep=""))
  par(mfrow=c(1,3),mar=c(4,4,4,4))
  
  if (length(diff_vals1)>4){
    dens <- apply(diff_vals1, 2, density, bw=2, na.rm=TRUE)
    plot(NA, xlim=c(-21,21), ylim=range(sapply(dens, "[", "y")),
      main=expression('1 ' < 'AGDD Z Score'),
      xlab='Model Bias (days)')
    mapply(lines, dens, col=1:length(dens),lwd=4)
    abline(v=0,lty=2)
  }
  if (length(diff_vals2)>4){
    dens <- apply(diff_vals2, 2, density, bw=2, na.rm=TRUE)
    plot(NA, xlim=c(-21,21), ylim=range(sapply(dens, "[", "y")),
      main=expression('1 ' <= 'AGDD Z Score < 2'),
      xlab='Model Bias (days)')
    mapply(lines, dens, col=1:length(dens),lwd=4)
    abline(v=0,lty=2)
  }
  if (length(diff_vals3)>4){
    dens <- apply(diff_vals3, 2, density, bw=2, na.rm=TRUE)
    plot(NA, xlim=c(-21,21), ylim=range(sapply(dens, "[", "y")),
      main=expression('AGDD Z Score > 2'),
      xlab='Model Bias (days)')
    mapply(lines, dens, col=1:length(dens),lwd=4)
    abline(v=0,lty=2)
    legend("topleft", legend=names(dens), fill=1:length(dens),bty='n')
  }
  
  dev.off()
}

