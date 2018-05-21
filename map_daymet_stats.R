library(ncdf4)
require(rgdal)
library(raster)

usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
usa_shp_proj <- spTransform(usa_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_crop <- crop(usa_shp_proj,extent(194214,2426712,-1660235,946090))

my.colors1 = colorRampPalette(c("white","yellow","orange","firebrick4"))
my.colors2 = colorRampPalette(c("white","red"))
my.colors3 = colorRampPalette(c("green4","green","white","deepskyblue","deepskyblue4"))
my.colors4 = colorRampPalette(c("firebrick4","orangered","orange","white",
  "cyan","deepskyblue","deepskyblue4"))
my.colors5 = colorRampPalette(c("red","white","blue"))

models <- c('Chill','Photo','SW_Jan01','SW_Mar13')
colors <- c('royalblue','springgreen','hotpink','orange')

sidelaps <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sidelap_mask.tif')
sidelaps_vals <- getValues(sidelaps)
w <- which(sidelaps_vals==0)

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("nobs*tif"),full.names=T,include.dirs=T,recursive=TRUE)
nobs <- raster(in_dirs_tile[1])
nobs[sidelaps==0] <- NA
pdf('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/nobs.pdf')
plot(nobs,main=c('Num. of Obs.'),
  col=my.colors2(15),colNA='black',axes=FALSE,zlim=c(20,34))
dev.off()

# Number of months < 5degC
below5 <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/below5_sum.tif')
below5_vals <- getValues(below5)
nobs_vals <- getValues(nobs)
below5_vals[is.na(nobs_vals)==1] <- NA
below5 <- setValues(below5,below5_vals)
plot(below5,main='Num. of Months T < 5C',col=my.colors5(9),colNA='black',axes=FALSE)

# RMSE
setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("RMSE*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile[])
s_vals <- getValues(s)
s_vals[s_vals>21] <- 21
s_vals[w,] <- NA
s <- setValues(s,s_vals)
pdf('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/rmse.pdf')
plot(s,main=c('Chill','Photo','SW-Jan1','SW-Mar13'),
  col=my.colors3(150),colNA='black',axes=FALSE)
dev.off()

rmse.group <- matrix(NA,9,4)
rmse.error.group <- matrix(NA,9,4)
for (i in 1:4){
  rmse.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=mean,na.rm=TRUE)[,2]
  rmse.error.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=sd,na.rm=TRUE)[,2]
}
x11(h=5,w=8)
#pdf(h=5,w=8,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_4/rmse_by_Tavg.pdf')
barCenters <- barplot(t(rmse.group),beside=TRUE,names=c(0,1,2,3,4,5,6,7,8),ylim=c(0,60),
  col=colors, ylab='RMSE (days)',xlab=expression('Months with T'[avg]*' < 5'~degree*'C'))
segments(barCenters, t(rmse.group)-t(rmse.error.group), barCenters,
  t(rmse.group)+t(rmse.error.group), lwd=1.5)
legend('topright',legend=models,
  cex=0.8,fill=colors,
  ncol=2)

# Correlation Coefficient
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("cor*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals[s_vals<0] <- 0 
s_vals[w,] <- NA
s <- setValues(s,s_vals)
pdf('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/cor.pdf')
plot(s,main=c('Chill','Photo','SW-Jan1','SW-Mar13'),
  col=my.colors4(150),colNA='black',zlim=c(0,1),axes=FALSE)
dev.off()

rmse.group <- matrix(NA,9,4)
rmse.error.group <- matrix(NA,9,4)
for (i in 1:4){
  rmse.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=mean,na.rm=TRUE)[,2]
  rmse.error.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=sd,na.rm=TRUE)[,2]
}
x11(h=5,w=8)
#pdf(h=5,w=8,'/projectnb/modislc/projects/te_phenology/figures/Storyboard_4/rmse_by_Tavg.pdf')
barCenters <- barplot(t(rmse.group),beside=TRUE,names=c(0,1,2,3,4,5,6,7,8),ylim=c(0,1),
  col=colors, ylab='RMSE (days)',xlab=expression('Months with T'[avg]*' < 5'~degree*'C'))
segments(barCenters, t(rmse.group)-t(rmse.error.group), barCenters,
  t(rmse.group)+t(rmse.error.group), lwd=1.5)
legend('topright',legend=models,
  cex=0.8,fill=colors,
  ncol=2)

# RMA Slope
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("slope*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals[s_vals<0] <- 0 
s_vals[s_vals>2] <- 2
s_vals[w,] <- NA
s <- setValues(s,s_vals)
plot(s,main=c('Chill','Photo','SW-Jan1','SW-Mar13'),
  col=my.colors3(150),colNA='black',zlim=c(0,2),axes=FALSE)

rmse.group <- matrix(NA,9,4)
rmse.error.group <- matrix(NA,9,4)
for (i in 1:4){
  rmse.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=mean,na.rm=TRUE)[,2]
  rmse.error.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=sd,na.rm=TRUE)[,2]
}
#11(h=5,w=8)
pdf(h=5,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/bar_cor_AGDDgt2sd.pdf')
barCenters <- barplot(t(rmse.group),beside=TRUE,names=c(0,1,2,3,4,5,6,7,8),ylim=c(0,1),
  col=colors, ylab='Correlation Coefficient',xlab=expression('Months with T'[avg]*' < 5'~degree*'C'))
segments(barCenters, t(rmse.group)-t(rmse.error.group), barCenters,
  t(rmse.group)+t(rmse.error.group), lwd=1.5)
legend('topright',legend=models,
  cex=0.8,fill=colors,
  ncol=4)
dev.off()

# Mean Bias Error
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("mbe*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals[s_vals < -21] <- -21 
s_vals[s_vals > 21] <- 21
s_vals[w,] <- NA
s <- setValues(s,s_vals)
plot(s,main=c('Chill','Photo','SW-Jan1','SW-Mar13'),
  col=my.colors4(150),colNA='black',zlim=c(-21,21),axes=FALSE)

# RMSE for pixel-year AGDD > 2SD
setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("rmse_2sd*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals[s_vals>21] <- 21
s_vals[w,] <- NA
s <- setValues(s,s_vals)
pdf('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/rmse_AGDDgt2sd.pdf')
plot(s,main=c('Chill','Photo','SW-Jan1','SW-Mar13'),
  col=my.colors3(150),colNA='black',axes=FALSE)
dev.off()

rmse.group <- matrix(NA,9,4)
rmse.error.group <- matrix(NA,9,4)
for (i in 1:4){
  rmse.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=mean,na.rm=TRUE)[,2]
  rmse.error.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=sd,na.rm=TRUE)[,2]
}
#x11(h=5,w=8)
pdf(h=5,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/bar_rmse_AGDDgt2sd.pdf')
barCenters <- barplot(t(rmse.group),beside=TRUE,names=c(0,1,2,3,4,5,6,7,8),ylim=c(0,27),
  col=colors, ylab='RMSE (days)',xlab=expression('Months with T'[avg]*' < 5'~degree*'C'))
segments(barCenters, t(rmse.group)-t(rmse.error.group), barCenters,
  t(rmse.group)+t(rmse.error.group), lwd=1.5)
legend('topright',legend=models,
  cex=0.8,fill=colors,
  ncol=4)
dev.off()