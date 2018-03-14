library(ncdf4)
require(rgdal)
library(raster)
library(foreach)
library(iterators)
library(doParallel)

source('~/Code/daymetr/R/download.daymet.tiles.r')

#args = commandArgs(trailingOnly=T) 
#tile_name = args[1]  

tile_name <- 'h17v02'

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',sep=''))
dir.create(paste(getwd(),'/CLIM/',sep=''))
dir.create(paste(getwd(),'/SHP/',sep=''))

H <- as.numeric(substring(tile_name,2,3))
V <- as.numeric(substring(tile_name,5,6))

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',sep=''))
ard_tiles <- readOGR('/projectnb/modislc/projects/landsat_sentinel/ARD/CONUS_ARD_grid/',
  'conus_ard_grid')
tile <- ard_tiles[which(ard_tiles$h==H & ard_tiles$v==V),]
# Determine lat/lon extent of tile for NED download
tile_latlon <- spTransform(tile,CRS("+proj=longlat +datum=WGS84"))
extent(tile_latlon)
tile_proj <- spTransform(tile,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
writeOGR(tile_proj,paste(getwd(),'/','SHP',sep=''),
  tile_name,driver="ESRI Shapefile",overwrite=TRUE)

tile <- readOGR(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/SHP/',sep=''),tile_name)
tile_proj <- spTransform(tile,CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')) 

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/CLIM/',sep=''))
download.daymet.tiles(lat1=tile_proj@bbox[2,2],lon1=tile_proj@bbox[1,1],
  lat2=tile_proj@bbox[2,1],lon2=tile_proj@bbox[1,2],start_yr=1981,end_yr=1981,param="tmax")

setwd('/projectnb/modislc/data/daymet/')
download.daymet.tiles(lat1=tile_proj@bbox[2,2],lon1=tile_proj@bbox[1,1],
  lat2=tile_proj@bbox[2,1],lon2=tile_proj@bbox[1,2],start_yr=1981,end_yr=2017,param="ALL")