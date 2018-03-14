require(raster)
require(rgdal)
require(gdalUtils)
require(rgeos)
require(spatialEco)
library(landsat)

# tile_name <- 'h23v14'
# H <- as.numeric(substring(tile_name,2,3))
# V <- as.numeric(substring(tile_name,5,6))
# 
# setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',sep=''))
# 
# ## GENERATE ARD TILE AND SAVE TO DIRECTORY
# ard_tiles <- readOGR('/projectnb/modislc/projects/landsat_sentinel/ARD/CONUS_ARD_grid/',
#   'conus_ard_grid')
# tile <- ard_tiles[which(ard_tiles$h==H & ard_tiles$v==V),]
# tile_proj <- spTransform(tile,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# 
# src_data <- '/projectnb/modislc/data/dem/usgs_ned/mosaic_dem.tif'
# DEM_proj <- gdalwarp(src_data,
#   dstfile='/projectnb/modislc/data/dem/usgs_ned/mosaic_dem_aea_proj.tif',
#   t_srs=projection(tile_proj),tr=c(30,30),r='bilinear',
#   output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)

src_data <- '/projectnb/modislc/data/dem/usgs_ned/mosaic_dem_aea_proj.tif'
DEM <- raster(src_data)

slope <- terrain(DEM,opt='slope',unit='degrees')
writeRaster(slope,filename='/projectnb/modislc/data/dem/usgs_ned/mosaic_slope_aea_proj.tif',
  format='GTiff',overwrite=TRUE)

aspect <- terrain(DEM,opt='aspect',unit='degrees')
writeRaster(aspect,filename='/projectnb/modislc/data/dem/usgs_ned/mosaic_aspect_aea_proj.tif',
  format='GTiff',overwrite=TRUE)
