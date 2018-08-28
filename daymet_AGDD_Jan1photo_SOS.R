##This script will generate long-term mean and annual 32km maps of spring and autumn phenology
##for a given Landast overlap scene

library(ncdf4)
require(rgdal)
library(raster)

args = commandArgs(trailingOnly=T)
tile_name = args[1]
#tile_name <- "h22v10"

#Load observed, 1km data
setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
load(file = "landsat2daymet")

#Native Daymet projection
daymet_crs <- CRS("+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +a=6378137 +rf=298.257223563 +lat_1=25 +lat_2=60")

#Import overlap shapefile
o <- readOGR(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
  tile_name,'/SHP/',sep=""),tile_name)

#Reproject overlap shapefile
o_reproj <- spTransform(o,daymet_crs)

AGDD <- matrix(NA,nrow(obs.SPR),36)
for (yr in 1982:2017){
  print(yr)
  
  setwd('/projectnb/modislc/data/climate/daymet')
  tmax <- brick(paste('daymet_v3_tmax_',yr,'_na.nc4',sep=''),var='tmax')
  tmin <- brick(paste('daymet_v3_tmin_',yr,'_na.nc4',sep=''),var='tmin')
  dayl <- brick(paste('daymet_v3_dayl_',yr,'_na.nc4',sep=''),var='dayl')

  tmax.crop <- crop(tmax,extent(o_reproj))
  tmin.crop <- crop(tmin,extent(o_reproj))
  dayl.crop <- crop(dayl,extent(o_reproj))
  
  tmean <- getValues((tmax.crop+tmin.crop)/2)
  tmean[tmean<0] <- 0
  
  photo <- getValues(dayl.crop)/3660
  tmean <- tmean*(photo/10)^6
  
  col <- (yr-1981)+3
  
  for (j in 1:nrow(obs.SPR)){
    if (is.na(obs.SPR[j,col])==0 & obs.SPR[j,col]>60){
      AGDD[j,(yr-1981)] <- sum(tmean[j,1:obs.SPR[j,col]])
    }
  }
}

assign('all.AGDD',AGDD)

#Calculate statistics for anomalously early observations
o_mean <- apply(all.AGDD,1,mean,na.rm=TRUE)
o_sd <- apply(all.AGDD,1,sd,na.rm=TRUE)
o_mean_all <- replicate(36, o_mean)
o_sd_all <- replicate(36, o_sd)
o_Zscore <- (all.AGDD-o_mean_all)/o_sd_all

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
save(o_mean,all.AGDD,o_Zscore,file = "daymet_AGDD_Jan1photo_exp6_SOS")