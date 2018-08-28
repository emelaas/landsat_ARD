library(ncdf4)
require(rgdal)
library(raster)

can_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/Canada/','Canada')
can_shp_proj <- spTransform(can_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
usa_shp_proj <- spTransform(usa_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_crop <- crop(usa_shp_proj,extent(194214,2426712,-1660235,946090))
can_crop <- crop(can_shp_proj,extent(194214,2426712,-1660235,946090))

my.colors1 = colorRampPalette(c("white","yellow","orange","firebrick4"))
my.colors2 = colorRampPalette(c("white","red"))
my.colors3 = colorRampPalette(c("green4","green","white","mediumpurple2","mediumorchid4"))
my.colors4 = colorRampPalette(c("firebrick4","orangered","orange","white",
  "cyan","deepskyblue","deepskyblue4"))
my.colors5 = colorRampPalette(c("red","white","blue"))

models <- c('Chill','Photo','SW Mar 17','SW Jan 1 (Mean AGDD)')
colors <- c('royalblue','springgreen','hotpink','orange')

sidelaps <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sidelap_mask.tif')
sidelaps_vals <- getValues(sidelaps)
w <- which(sidelaps_vals==0)

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("nobs*tif"),full.names=T,include.dirs=T,recursive=TRUE)
nobs <- raster(in_dirs_tile[1])
nobs[sidelaps==0] <- NA
pdf('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/nobs.pdf')
plot(usa_crop,col='white',main='Number of Annual Phenology Obs.')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray30',add=TRUE)
plot(can_crop,col='gray30',add=TRUE)
plot(nobs,main=c('Num. of Obs.'),
  col=my.colors2(15),axes=FALSE,zlim=c(20,34),add=TRUE)
dev.off()

# Number of months < 5degC
below5 <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/below5_sum.tif')
below5_vals <- getValues(below5)
nobs_vals <- getValues(nobs)
below5_vals[is.na(nobs_vals)==1] <- NA
below5 <- setValues(below5,below5_vals)
pdf('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/below5_sum.pdf')
plot(usa_crop,col='white',main='Num. of Months T < 5C')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
plot(below5,col=my.colors5(9),axes=FALSE,add=TRUE)
dev.off()

# RMSE
setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("RMSE*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile[])
s_vals <- getValues(s)
s_vals[s_vals>14] <- 14
s_vals[w,] <- NA
s <- setValues(s,s_vals)
#x11(h=7,w=7)
pdf(h=7,w=7,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/rmse.pdf')
par(mfrow=c(2,2),mar=c(2,0,1,0))
for (i in 1:4){
  plot(usa_crop,col='white',main=models[i])
  rect(194214,-1660235,2426712,946090,col='lightblue')
  plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
  if (i==1) plot(s[[i]],col=my.colors3(150),axes=FALSE,add=TRUE)
  if (i!=1) plot(s[[i]],col=my.colors3(150),axes=FALSE,add=TRUE,legend=FALSE)
}
dev.off()

rmse.group <- matrix(NA,9,4)
rmse.error.group <- matrix(NA,9,4)
for (i in 1:4){
  rmse.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=mean,na.rm=TRUE)[,2]
  rmse.error.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=sd,na.rm=TRUE)[,2]
}
#x11(h=5,w=8)
pdf(h=5,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/rmse_by_Tavg.pdf')
barCenters <- barplot(t(rmse.group),beside=TRUE,names=c(0,1,2,3,4,5,6,7,8),ylim=c(0,14),
  col=colors, ylab='RMSE (days)',xlab=expression('Months with T'[avg]*' < 5'~degree*'C'))
segments(barCenters, t(rmse.group)-t(rmse.error.group), barCenters,
  t(rmse.group)+t(rmse.error.group), lwd=1.5)
legend('topleft',legend=models,
  cex=0.8,fill=colors,
  ncol=1)
dev.off()

# Correlation Coefficient
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("cor*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals[s_vals<0] <- 0 
s_vals[w,] <- NA
s <- setValues(s,s_vals)
#x11(h=7,w=7)
pdf(h=7,w=7,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/cor.pdf')
par(mfrow=c(2,2),mar=c(2,0,1,0))
for (i in 1:4){
  plot(usa_crop,col='white',main=models[i])
  rect(194214,-1660235,2426712,946090,col='lightblue')
  plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
  if (i==1) plot(s[[i]],col=my.colors3(150),axes=FALSE,add=TRUE)
  if (i!=1) plot(s[[i]],col=my.colors3(150),axes=FALSE,add=TRUE,legend=FALSE)
}
dev.off()

cor.group <- matrix(NA,9,4)
cor.error.group <- matrix(NA,9,4)
for (i in 1:4){
  cor.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=mean,na.rm=TRUE)[,2]
  cor.error.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=sd,na.rm=TRUE)[,2]
}
#x11(h=5,w=8)
pdf(h=5,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/cor_by_Tavg.pdf')
barCenters <- barplot(t(cor.group),beside=TRUE,names=c(0,1,2,3,4,5,6,7,8),ylim=c(0,1),
  col=colors, ylab='Correlation Coefficient',xlab=expression('Months with T'[avg]*' < 5'~degree*'C'))
segments(barCenters, t(cor.group)-t(cor.error.group), barCenters,
  t(cor.group)+t(cor.error.group), lwd=1.5)
legend('topright',legend=models,
  cex=0.8,fill=colors,
  ncol=3)
dev.off()

# RMA Slope
slope.group <- matrix(NA,9,4)
slope.error.group <- matrix(NA,9,4)
for (i in 1:4){
  slope.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=mean,na.rm=TRUE)[,2]
  slope.error.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=sd,na.rm=TRUE)[,2]
}
#x11(h=5,w=8)
pdf(h=5,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/slope_by_Tavg.pdf')
barCenters <- barplot(t(slope.group),beside=TRUE,names=c(0,1,2,3,4,5,6,7,8),ylim=c(0,2),
  col=colors, ylab='RMA Slope',xlab=expression('Months with T'[avg]*' < 5'~degree*'C'))
segments(barCenters, t(slope.group)-t(slope.error.group), barCenters,
  t(slope.group)+t(slope.error.group), lwd=1.5)
legend('topright',legend=models,
  cex=0.8,fill=colors,
  ncol=1)
dev.off()

# Mean Bias Error
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("mbe*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals[s_vals < -14] <- -14
s_vals[s_vals > 14] <- 14
s_vals[w,] <- NA
s <- setValues(s,s_vals)
#x11(h=7,w=7)
pdf(h=7,w=7,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/mbe.pdf')
par(mfrow=c(2,2),mar=c(2,0,1,0))
for (i in 1:4){
  plot(usa_crop,col='white',main=models[i])
  rect(194214,-1660235,2426712,946090,col='lightblue')
  plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
  if (i==1) plot(s[[i]],col=my.colors3(150),axes=FALSE,add=TRUE)
  if (i!=1) plot(s[[i]],col=my.colors3(150),axes=FALSE,add=TRUE,legend=FALSE)
}
dev.off()

# RMSE for pixel-year AGDD > 2SD
setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("rmse_2sd*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals[s_vals>14] <- 14
s_vals[w,] <- NA
s <- setValues(s,s_vals)
#x11(h=7,w=7)
pdf(h=7,w=7,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/rmse_AGDDgt2sd.pdf')
par(mfrow=c(2,2),mar=c(2,0,1,0))
for (i in 1:4){
  plot(usa_crop,col='white',main=models[i])
  rect(194214,-1660235,2426712,946090,col='lightblue')
  plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
  if (i==1) plot(s[[i]],col=my.colors3(150),axes=FALSE,add=TRUE)
  if (i!=1) plot(s[[i]],col=my.colors3(150),axes=FALSE,add=TRUE,legend=FALSE)
}
dev.off()

rmse.group <- matrix(NA,9,4)
rmse.error.group <- matrix(NA,9,4)
for (i in 1:4){
  rmse.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=mean,na.rm=TRUE)[,2]
  rmse.error.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=sd,na.rm=TRUE)[,2]
}
#x11(h=5,w=8)
pdf(h=5,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/bar_rmse_AGDDgt2sd.pdf')
barCenters <- barplot(t(rmse.group),beside=TRUE,names=c(0,1,2,3,4,5,6,7,8),ylim=c(0,27),
  col=colors, ylab='RMSE (days)',xlab=expression('Months with T'[avg]*' < 5'~degree*'C'))
segments(barCenters, t(rmse.group)-t(rmse.error.group), barCenters,
  t(rmse.group)+t(rmse.error.group), lwd=1.5)
legend('topright',legend=models,
  cex=0.8,fill=colors,
  ncol=1)
dev.off()


# Load phenocam transition dates and associated Daymet time series
#load("~/Code/GitHub/phenor/data/phenocam_DB.rda") #old version from KH phenor package
load("/projectnb/modislc/users/emelaas/MSB/PhenoCam/phenocam_DB.RData")

# concat locations data into a matrix with the first row
# being the latitude and the second longitude
location = do.call("cbind",lapply(phenocam_DB,function(x){
  if(!is.null(x)){
    matrix(rep(x$location, ncol(x$Ti)), 2, ncol(x$Ti))
  }
}))
loc <- data.frame(location[2,],location[1,])
colnames(loc) <- c('Lon','Lat')
coordinates(loc) <- ~Lon+Lat
proj4string(loc) <- CRS("+proj=longlat +datum=WGS84")
loc_proj <- spTransform(loc, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))

pdf(h=8,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/pcam_locs.pdf')
plot(usa_crop,col='white',main='DBF PhenoCam Sites')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
plot(loc_proj,add=TRUE,pch=21,bg='yellow')
dev.off()