require(raster)
require(rgdal)
require(gdalUtils)
require(rgeos)
require(spatialEco)
library(landsat)

## GENERATE DEM MOSAIC AND CLIP TO ARD TILE
# Specify a for loop to create a list object, containing raster objects
setwd('/projectnb/modislc/data/dem/usgs_ned')
rasters1 <- list.files(path=getwd(),pattern=glob2rx("*img"),
  full.names=T,include.dirs=T,recursive=T)
rast.list <- list()
for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }

# Create mosaic using do.call on the list of raster objects
rast.list$fun <- mean
rast.mosaic <- do.call(mosaic,rast.list)
writeRaster(rast.mosaic,filename=paste(getwd(),'/mosaic',sep=""),
  format='GTiff',overwrite=TRUE)



