## This script generates EVI2 GeoTIFF for each image in the ARD stack
## Still working on topographic correction component (need to modify for each LC)

require(raster)
require(rgdal)
require(gdalUtils)
require(rgeos)
require(jsonlite)
require(spatialEco)
library(foreach)
library(iterators)
library(doParallel)
library(landsat)

#Register the parallel backend
registerDoParallel(16)

source(file='/usr3/graduate/emelaas/Code/R/landsat_sentinel/v1_3/topocorr_v2.R')

# args = commandArgs(trailingOnly=T) 
# tile_name = args[1]  


tile_name <- 'h17v15'
H <- as.numeric(substring(tile_name,2,3))
V <- as.numeric(substring(tile_name,5,6))

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',sep=''))

## GENERATE ARD TILE AND SAVE TO DIRECTORY
ard_tiles <- readOGR('/projectnb/modislc/projects/landsat_sentinel/ARD/CONUS_ARD_grid/',
  'conus_ard_grid')
tile <- ard_tiles[which(ard_tiles$h==H & ard_tiles$v==V),]
# Determine lat/lon extent of tile for NED download
tile_latlon <- spTransform(tile,CRS("+proj=longlat +datum=WGS84"))
extent(tile_latlon)
tile_proj <- spTransform(tile,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
writeOGR(tile_proj,paste(getwd(),'/','SHP',sep=''),
  tile_name,driver="ESRI Shapefile",overwrite=TRUE)

# ## GENERATE DEM MOSAIC AND CLIP TO ARD TILE
# # Specify a for loop to create a list object, containing raster objects
# setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/DEM',sep=''))
# rasters1 <- list.files(path=getwd(),pattern=glob2rx("*img"),
#   full.names=T,include.dirs=T,recursive=T)
# rast.list <- list()
# for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }
# 
# # Create mosaic using do.call on the list of raster objects
# rast.list$fun <- mean
# rast.mosaic <- do.call(mosaic,rast.list)
# writeRaster(rast.mosaic,filename=paste(getwd(),'/mosaic',sep=""),
#   format='GTiff',overwrite=TRUE)
# 
# setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/SHP',sep=''))
# DEM <- gdalwarp(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/DEM/mosaic.tif',sep=''),
#   dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/DEM/mosaic_tile.tif',sep=''),
#   t_srs=projection(tile_proj),ts=c(5000,5000),
#   cutline=paste(tile_name,'.shp',sep=''),cl=tile_name,
#   crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/SHP',sep=''))
src_data2 <- '/projectnb/modislc/data/dem/usgs_ned/mosaic.tif'
DEM <- gdalwarp(src_data2,
  dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/DEM/mosaic_tile.tif',sep=''),
  t_srs=projection(tile_proj),ts=c(5000,5000),
  cutline=paste(tile_name,'.shp',sep=''),cl=tile_name,
  crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)

slope <- terrain(DEM,opt='slope',unit='degrees')
aspect <- terrain(DEM,opt='aspect',unit='degrees')

src_data2 <- '/projectnb/modislc/data/lc_database/regional/united_states/NLCD2006_landcover_4-20-11_se5/nlcd2006_landcover_4-20-11_se5.img'
LC <- gdalwarp(src_data2,dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/NLCD/nlcd.tif',sep=''),
  t_srs=projection(tile_proj),ts=c(5000,5000),
  cutline=paste(tile_name,'.shp',sep=''),cl=tile_name,
  crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)


# LOAD IN SURFACE REFLECTANCE BANDS, QA LAYER, SOLAR ZENITH & AZIMUTH

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/IMG',sep=''))
in_dirs <- list.files(path=getwd(),pattern=glob2rx("L*"),
  full.names=T,include.dirs=T)

all_count <- foreach(i = 1:length(in_dirs), .combine = rbind) %dopar% {
  print(i)

  # Landsat 8
  if (substr(in_dirs[i],64,64)==8){
    nir <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB5.tif',sep='')) # near infrared reflectance
    red <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB4.tif',sep='')) # red reflectance
    QA <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_PIXELQA.tif',sep='')) # pixel QA

    s <- stack(red,nir,QA)
    band_vals <- getValues(s)
    band_vals[band_vals<=0] <- NA

    qa <- band_vals[,3]
    w <- which(qa!=322 & qa!=386 & qa!=834 & qa!=898 & qa!=1346)
    band_vals[w,1:2] <- NA

    # Landsat 4-7
  } else {
    nir <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB4.tif',sep='')) # near infrared reflectance
    red <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB3.tif',sep='')) # red reflectance
    QA <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_PIXELQA.tif',sep='')) # pixel QA

    s <- stack(red,nir,QA)
    band_vals <- getValues(s)
    band_vals[band_vals<=0] <- NA

    qa <- band_vals[,3]
    w <- which(qa!=66 & qa!=130)
    band_vals[w,1:2] <- NA
  }

#   # TOPOGRAPHIC CORRECTION (NEED TO MODIFY TO ADJUST FOR DIFFERENT LC TYPES)
#   SAA <- 0.01*raster(paste(in_dirs[i],'/',substr(in_dirs[i],57,88),
#     '_C01_V01_SOA4.tif',sep='')) # solar azimuth angle
#   SZA <- 0.01*raster(paste(in_dirs[i],'/',substr(in_dirs[i],57,88),
#     '_C01_V01_SOZ4.tif',sep='')) # solar zenith angle
#
#   if (length(which(is.na(band_vals[,1])==0)) > 4){
#     b3_corr <- as.numeric(t(topocorr_v2(x=as.matrix(setValues(b3,band_vals[,1])),
#       slope=as.matrix(slope), aspect=as.matrix(aspect), sunelev=90-as.matrix(SZA),
#       sunazimuth=as.matrix(SAA), method='rotational')))
#     b4_corr <- as.numeric(t(topocorr_v2(x=as.matrix(setValues(b4,band_vals[,2])),
#       slope=as.matrix(slope), aspect=as.matrix(aspect), sunelev=90-as.matrix(SZA),
#       sunazimuth=as.matrix(SAA), method='rotational')))
#   } else {
#     b3_corr <- as.numeric(matrix(NA,ncell(slope),1))
#     b4_corr <- as.numeric(matrix(NA,ncell(slope),1))
#   }

  evi2 <- 2.5*(band_vals[,2]/10000 - band_vals[,1]/10000)/(band_vals[,2]/10000 + 2.4*band_vals[,1]/10000 + 1)
  evi2_map <- setValues(nir,round(evi2*10000))
  writeRaster(evi2_map,filename=paste(in_dirs[i],'/evi2',sep=""),
    format='GTiff',overwrite=TRUE)

  count <- i
}


