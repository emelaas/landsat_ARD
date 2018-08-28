library(ncdf4)
require(rgdal)
library(raster)

args = commandArgs(trailingOnly=T)
year = as.numeric(args[1])
# m = as.numeric(args[2])

# year = 1998
m = 1

print(year)

# Find all Daymet tiles overlapping with ARD tile
setwd('/projectnb/modislc/data/climate/daymet')
tmax <- raster('daymet_v3_tmax_1981_na.nc4',varname='tmax')
cells <- setValues(tmax,seq(1,ncell(tmax)))

usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
usa_shp_proj <- spTransform(usa_shp,
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_crop <- crop(usa_shp_proj,extent(194214,2426712,-1660235,946090))

tiles <- read.table('~/Code/GitHub/ard_tiles.txt')


tmp1  <- matrix(NA,ncell(cells),1)
tmp2 <- matrix(NA,ncell(cells),1)
tmp3 <- matrix(NA,ncell(cells),1)
tmp4 <- matrix(NA,ncell(cells),1)
tmp5 <- matrix(NA,ncell(cells),1)

for (i in 1:170){
  print(c(year,i))
  
  tile_name <- as.character(tiles[i,1])
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/PHENO_1KM/',sep=''))
  load(file = "daymet_AGDD_Jan1"); rm('all.AGDD','o_Zscore')
  load(file = "daymet_AGDD")
  load(file = paste("daymet_predict_v5_",m,sep=""))
  load(file = "landsat2daymet")
  load(file = "landsat2daymet_aut")
  all.oSPR <- obs.SPR[,4:39]
  all.oAUT <- obs.AUT[,4:39]
  
  print(nrow(all.pSPR))
  
  w <- which(is.na(tmp1)==0)
  w2 <- which(obs.SPR[,1] %in% w)
  all.pSPR[w2,(year-1981)] <- tmp1[obs.SPR[w2,1]]
  tmp1[obs.SPR[,1],1] <- all.pSPR[,(year-1981)]
  
  if (m==1){
    w <- which(is.na(tmp2)==0)
    w2 <- which(obs.SPR[,1] %in% w)
    o_Zscore[w2,(year-1981)] <- tmp2[obs.SPR[w2,1]]
    tmp2[obs.SPR[,1],1] <- o_Zscore[,(year-1981)]
    
    # Observed Spring
    w <- which(is.na(tmp3)==0)
    w2 <- which(obs.SPR[,1] %in% w)
    all.oSPR[w2,(year-1981)] <- tmp3[obs.SPR[w2,1]]
    tmp3[obs.SPR[,1],1] <- all.oSPR[,(year-1981)]
    
    # Observed Autumn
    w <- which(is.na(tmp4)==0)
    w2 <- which(obs.SPR[,1] %in% w)
    all.oAUT[w2,(year-1981)] <- tmp4[obs.SPR[w2,1]]
    tmp4[obs.SPR[,1],1] <- all.oAUT[,(year-1981)]

    # Mean AGDD between Jan 1 and mean SOS
    w <- which(is.na(tmp5)==0)
    w2 <- which(obs.SPR[,1] %in% w)
    o_mean[w2] <- tmp5[obs.SPR[w2,1]]
    tmp5[obs.SPR[,1],1] <- o_mean
  }
}

sidelaps <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sidelap_mask.tif')
sprLTM <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sprLTM.tif')

# spr_mod <- setValues(cells,tmp1)
# spr_mod <- crop(spr_mod,usa_crop)
# spr_mod[sidelaps==0] <- NA
# spr_mod[is.na(sprLTM)==1] <- NA
# writeRaster(spr_mod,
#   filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr_mod',m,"_",year,sep=""),
#   format='GTiff',overwrite=TRUE)

if (m == 1){
#   Zscore <- setValues(cells,tmp2)
#   Zscore <- crop(Zscore,usa_crop)
#   Zscore[sidelaps==0] <- NA
#   Zscore[is.na(sprLTM)==1] <- NA
#   writeRaster(Zscore,
#     filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/agdd_Zscore_',year,sep=""),
#     format='GTiff',overwrite=TRUE)
  
  spr_obs <- setValues(cells,tmp3)
  spr_obs <- crop(spr_obs,usa_crop)
  #spr_obs[sidelaps==0] <- NA
  #spr_obs[is.na(sprLTM)==1] <- NA
  writeRaster(spr_obs,
    filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr_wsidelaps',year,sep=""),
    format='GTiff',overwrite=TRUE)
  
  aut_obs <- setValues(cells,tmp4)
  aut_obs <- crop(aut_obs,usa_crop)
  #aut_obs[sidelaps==0] <- NA
  #aut_obs[is.na(sprLTM)==1] <- NA
  writeRaster(aut_obs,
    filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/aut_wsidelaps',year,sep=""),
    format='GTiff',overwrite=TRUE)
  
#   agdd_Jan1 <- setValues(cells,tmp5)
#   agdd_Jan1 <- crop(agdd_Jan1,usa_crop)
#   agdd_Jan1[sidelaps==0] <- NA
#   agdd_Jan1[is.na(sprLTM)==1] <- NA
#   writeRaster(agdd_Jan1,
#     filename='/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/agdd_Jan1',
#     format='GTiff',overwrite=TRUE)
}


