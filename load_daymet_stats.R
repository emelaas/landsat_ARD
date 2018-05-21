library(ncdf4)
require(rgdal)
library(raster)

# args = commandArgs(trailingOnly=T)
# m = as.numeric(args[1])
m <- 1

usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
usa_shp_proj <- spTransform(usa_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_crop <- crop(usa_shp_proj,extent(194214,2426712,-1660235,946090))

# Find all Daymet tiles overlapping with ARD tile
setwd('/projectnb/modislc/data/daymet')
tmax <- raster('daymet_v3_tmax_1981_na.nc4',varname='tmax')
cells <- setValues(tmax,seq(1,ncell(tmax)))

tiles <- read.table('~/Code/GitHub/ard_tiles.txt')

rmse_map <- matrix(NA,ncell(cells),1)
mbe_map <- matrix(NA,ncell(cells),1)
cor_map <- matrix(NA,ncell(cells),1)
slope_map <- matrix(NA,ncell(cells),1)
rmse_2sd_map <- matrix(NA,ncell(cells),1)
nobs_map <- matrix(NA,ncell(cells),1)
below5_map <- matrix(NA,ncell(cells),1)
obsSD_map <- matrix(NA,ncell(cells),1)
agddSD_map <- matrix(NA,ncell(cells),1)
IQR_map <- matrix(NA,ncell(cells),1)
for (i in 1:172){

  print(c(i,m))
  
  tile_name <- tiles[i,1]
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=''))
  load(file = "landsat2daymet")
  load(file = "below5_sum")
  load(file = "daymet_AGDD")
  load(file = paste("daymet_predict_",m,sep=""))
  
  nobs.pix[nobs.pix==0] <- NA
  
  #o_sd <- apply(obs.SPR[,4:39],1,sd,na.rm=TRUE)
  agdd_sd <- apply(all.AGDD,1,sd,na.rm=TRUE)
  #IQR <- apply(obs.SPR[,4:39],1,quantile,0.75,1,na.rm=TRUE)-apply(obs.SPR[,4:39],1,quantile,0.25,1,na.rm=TRUE)
  
  w <- which(is.na(nobs_map)==0)
  w2 <- which(obs.SPR[,1] %in% w)
  
  rmse.pix[w2,1] <- rmse_map[obs.SPR[w2,1]]
  rmse_map[obs.SPR[,1]] <- rmse.pix
  
  mbe.pix[w2,1] <- mbe_map[obs.SPR[w2,1]]
  mbe_map[obs.SPR[,1]] <- mbe.pix
  
  cor.pix[w2,1] <- cor_map[obs.SPR[w2,1]]
  cor_map[obs.SPR[,1]] <- cor.pix[,1]
  
  slope.pix[w2,1] <- slope_map[obs.SPR[w2,1]]
  slope_map[obs.SPR[,1]] <- slope.pix
  
  rmse_2sd.pix[w2,1] <- rmse_2sd_map[obs.SPR[w2,1]]
  rmse_2sd_map[obs.SPR[,1]] <- rmse_2sd.pix
  
  nobs.pix[w2,1] <- nobs_map[obs.SPR[w2,1]]
  nobs_map[obs.SPR[,1]] <- nobs.pix
  
  below5_sum[w2] <- below5_map[obs.SPR[w2,1]]
  below5_map[obs.SPR[,1]] <- below5_sum
  
  #o_sd[w2] <- obsSD_map[obs.SPR[w2,1]]
  #obsSD_map[obs.SPR[,1]] <- o_sd
  
  agdd_sd[w2] <- agddSD_map[obs.SPR[w2,1]]
  agddSD_map[obs.SPR[,1]] <- agdd_sd
  
  #IQR[w2] <- IQR_map[obs.SPR[w2,1]]
  #IQR_map[obs.SPR[,1]] <- IQR
}

sidelaps <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sidelap_mask.tif')

rmse_mod <- setValues(cells,rmse_map)
rmse_mod <- crop(rmse_mod,usa_crop)
rmse_mod[sidelaps==0] <- NA
writeRaster(rmse_mod,
  filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/RMSE',m,sep=""),
  format='GTiff',overwrite=TRUE)

mbe_mod <- setValues(cells,mbe_map)
mbe_mod <- crop(mbe_mod,usa_crop)
mbe_mod[sidelaps==0] <- NA
writeRaster(mbe_mod,
  filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/mbe',m,sep=""),
  format='GTiff',overwrite=TRUE)

cor_mod <- setValues(cells,cor_map)
cor_mod <- crop(cor_mod,usa_crop)
cor_mod[sidelaps==0] <- NA
writeRaster(cor_mod,
  filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/cor',m,sep=""),
  format='GTiff',overwrite=TRUE)

slope_mod <- setValues(cells,slope_map)
slope_mod <- crop(slope_mod,usa_crop)
slope_mod[sidelaps==0] <- NA
writeRaster(slope_mod,
  filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/slope',m,sep=""),
  format='GTiff',overwrite=TRUE)

rmse_2sd_mod <- setValues(cells,rmse_2sd_map)
rmse_2sd_mod <- crop(rmse_2sd_mod,usa_crop)
rmse_2sd_mod[sidelaps==0] <- NA
writeRaster(rmse_2sd_mod,
  filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/rmse_2sd_',m,sep=""),
  format='GTiff',overwrite=TRUE)

nobs_mod <- setValues(cells,nobs_map)
nobs_mod <- crop(nobs_mod,usa_crop)
nobs_mod[sidelaps==0] <- NA
writeRaster(nobs_mod,
  filename='/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/nobs',
  format='GTiff',overwrite=TRUE)

below5_mod <- setValues(cells,below5_map)
below5_mod <- crop(below5_mod,usa_crop)
below5_mod[sidelaps==0] <- NA
writeRaster(below5_mod,
  filename='/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/below5_sum',
  format='GTiff',overwrite=TRUE)

obsSD_mod <- setValues(cells,obsSD_map)
obsSD_mod <- crop(obsSD_mod,usa_crop)
obsSD_mod[sidelaps==0] <- NA
writeRaster(obsSD_mod,
  filename='/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/obsSD',
  format='GTiff',overwrite=TRUE)

agddSD_mod <- setValues(cells,agddSD_map)
agddSD_mod <- crop(agddSD_mod,usa_crop)
agddSD_mod[sidelaps==0] <- NA
writeRaster(agddSD_mod,
  filename='/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/agddSD',
  format='GTiff',overwrite=TRUE)

# IQR_mod <- setValues(cells,IQR_map)
# IQR_mod <- crop(IQR_mod,usa_crop)
# IQR_mod[sidelaps==0] <- NA
# writeRaster(IQR_mod,
#   filename='/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/IQR_25_75',
#   format='GTiff',overwrite=TRUE)