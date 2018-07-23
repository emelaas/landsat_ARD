library(ncdf4)
require(rgdal)
library(raster)

args = commandArgs(trailingOnly=T)
year = as.numeric(args[1])

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

for (m in 1:4){

  tmp <- matrix(NA,ncell(cells),1)
  tmp2 <- matrix(NA,ncell(cells),1)
  tmp3 <- matrix(NA,ncell(cells),1)

  for (i in 1:172){
    print(c(year,i))

    tile_name <- as.character(tiles[i,1])

    setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
      tile_name,'/PHENO_1KM/',sep=''))
    load(file = "daymet_AGDD")
    load(file = paste("daymet_predict_v5_",m,sep=""))
    load(file = "landsat2daymet")
    all.oSPR <- obs.SPR[,4:39]

    w <- which(is.na(tmp)==0)
    w2 <- which(obs.SPR[,1] %in% w)
    all.pSPR[w2,(year-1981)] <- tmp[obs.SPR[w2,1]]
    tmp[obs.SPR[,1],1] <- all.pSPR[,(year-1981)]

    w <- which(is.na(tmp2)==0)
    w2 <- which(obs.SPR[,1] %in% w)
    o_Zscore[w2,(year-1981)] <- tmp2[obs.SPR[w2,1]]
    tmp2[obs.SPR[,1],1] <- o_Zscore[,(year-1981)]

    w <- which(is.na(tmp3)==0)
    w2 <- which(obs.SPR[,1] %in% w)
    all.oSPR[w2,(year-1981)] <- tmp3[obs.SPR[w2,1]]
    tmp3[obs.SPR[,1],1] <- all.oSPR[,(year-1981)]
  }

  sidelaps <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sidelap_mask.tif')
  sprLTM <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sprLTM.tif')

  spr_mod <- setValues(cells,tmp)
  spr_mod <- crop(spr_mod,usa_crop)
  spr_mod[sidelaps==0] <- NA
  spr_mod[is.na(sprLTM)==1] <- NA
  writeRaster(spr_mod,
    filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr_mod',m,"_",year,sep=""),
    format='GTiff',overwrite=TRUE)

  Zscore <- setValues(cells,tmp2)
  Zscore <- crop(Zscore,usa_crop)
  Zscore[sidelaps==0] <- NA
  Zscore[is.na(sprLTM)==1] <- NA
  writeRaster(Zscore,
    filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/agdd_Zscore_',year,sep=""),
    format='GTiff',overwrite=TRUE)

  spr_obs <- setValues(cells,tmp3)
  spr_obs <- crop(spr_obs,usa_crop)
  spr_obs[sidelaps==0] <- NA
  spr_obs[is.na(sprLTM)==1] <- NA
  writeRaster(spr_obs,
    filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/spr',year,sep=""),
    format='GTiff',overwrite=TRUE)
}
