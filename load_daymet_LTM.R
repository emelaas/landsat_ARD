library(ncdf4)
require(rgdal)
library(raster)

# Find all Daymet tiles overlapping with ARD tile
setwd('/projectnb/modislc/data/daymet')
tmax <- raster('daymet_v3_tmax_1981_na.nc4',varname='tmax')
cells <- setValues(tmax,seq(1,ncell(tmax)))

usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
usa_shp_proj <- spTransform(usa_shp,
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_crop <- crop(usa_shp_proj,extent(194214,2426712,-1660235,946090))

tiles <- read.table('~/Code/GitHub/ard_tiles.txt')

tmp <- matrix(NA,ncell(cells),1)

for (i in 1:172){
  print(i)

  tile_name <- as.character(tiles[i,1])

  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/PHENO_1KM/',sep=''))
  load(file = "landsat2daymet")

  w <- which(is.na(tmp)==0)
  w2 <- which(obs.SPR[,1] %in% w)
  obs.SPR[w2,2] <- tmp[obs.SPR[w2,1]]
  tmp[obs.SPR[,1],1] <- obs.SPR[,2]
}

sidelaps <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sidelap_mask.tif')

sprLTM <- setValues(cells,tmp)
sprLTM <- crop(sprLTM,usa_crop)
sprLTM[sidelaps==0] <- NA
writeRaster(sprLTM,
  filename='/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sprLTM',
  format='GTiff',overwrite=TRUE)
