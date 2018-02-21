library(ncdf4)
require(rgdal)
library(raster)
library(foreach)
library(iterators)
library(doParallel)

source('~/Code/daymetr/R/download.daymet.tiles.r')

args = commandArgs(trailingOnly=T) 
tile_name = args[1]  

#tile_name <- 'h22v15'

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',sep=''))
dir.create(paste(getwd(),'/CLIM/',sep=''))

tile <- readOGR(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/SHP/',sep=''),tile_name)
tile_proj <- spTransform(tile,CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')) 

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/CLIM/',sep=''))
download.daymet.tiles(lat1=tile_proj@bbox[2,2],lon1=tile_proj@bbox[1,1],
  lat2=tile_proj@bbox[2,1],lon2=tile_proj@bbox[1,2],start_yr=1981,end_yr=2016,param="ALL")