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

source(file='/usr3/graduate/emelaas/Code/R/landsat_sentinel/v1_3/topocorr_v3.R')

#args = commandArgs(trailingOnly=T) 
#tile_name = args[1]  

tile_name <- 'h29v07'

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

## LOAD IN SURFACE REFLECTANCE BANDS, QA LAYER, SOLAR ZENITH & AZIMUTH

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/IMG',sep=''))
in_dirs <- list.files(path=getwd(),pattern=glob2rx("L*"),
  full.names=T,include.dirs=T)

all_count <- foreach(i = 1:length(in_dirs), .combine = rbind) %dopar% {
  
  # Landsat 8
  if (substr(in_dirs[i],64,64)==8){
    blue <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB2.tif',sep='')) # blue reflectance
    nir <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB5.tif',sep='')) # near infrared reflectance
    red <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB4.tif',sep='')) # red reflectance
    swir1 <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB6.tif',sep='')) # shortwave infrared reflectance
    swir2 <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB7.tif',sep='')) # shortwave infrared reflectance
    QA <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_PIXELQA.tif',sep='')) # pixel QA
    
    s <- stack(blue,red,nir,swir1,swir2,QA)
    band_vals <- getValues(s)
    band_vals[band_vals<=0] <- NA
    
    qa <- band_vals[,6]
    w <- which(qa!=322 & qa!=386 & qa!=834 & qa!=898 & qa!=1346)
    band_vals[w,1:5] <- NA
    
    # Landsat 4-7
  } else {
    blue <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB1.tif',sep='')) # blue reflectance
    nir <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB4.tif',sep='')) # near infrared reflectance
    red <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB3.tif',sep='')) # red reflectance
    swir1 <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB5.tif',sep='')) # shortwave infrared reflectance
    swir2 <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_SRB7.tif',sep='')) # shortwave infrared reflectance
    QA <- raster(paste(in_dirs[i],'/',substr(in_dirs[i],61,92),
      '_C01_V01_PIXELQA.tif',sep='')) # pixel QA
    
    s <- stack(blue,red,nir,swir1,swir2,QA)
    band_vals <- getValues(s)
    band_vals[band_vals<=0] <- NA
    
    qa <- band_vals[,6]
    w <- which(qa!=66 & qa!=130)
    band_vals[w,1:5] <- NA
  }
  
#   evi2 <- 2.5*(band_vals[,2]/10000 - band_vals[,1]/10000)/(band_vals[,2]/10000 + 2.4*band_vals[,1]/10000 + 1)
#   evi2_map <- setValues(nir,round(evi2*10000))
#   writeRaster(evi2_map,filename=paste(in_dirs[i],'/evi2',sep=""),
#     format='GTiff',overwrite=TRUE)
#   
#   ndvi <- (band_vals[,2]/10000 - band_vals[,1]/10000)/(band_vals[,2]/10000 + band_vals[,1]/10000)
#   ndvi_map <- setValues(nir,round(ndvi*10000))
#   writeRaster(ndvi_map,filename=paste(in_dirs[i],'/ndvi',sep=""),
#     format='GTiff',overwrite=TRUE)
#   
#   ndwi <- (band_vals[,2]/10000 - band_vals[,4]/10000)/(band_vals[,2]/10000 + band_vals[,4]/10000)
#   ndwi_map <- setValues(nir,round(ndwi*10000))
#   writeRaster(ndwi_map,filename=paste(in_dirs[i],'/ndwi',sep=""),
#     format='GTiff',overwrite=TRUE)
#   
#   msi <- (band_vals[,3]/10000)/(band_vals[,2]/10000)
#   msi_map <- setValues(nir,round(msi*10000))
#   writeRaster(msi_map,filename=paste(in_dirs[i],'/msi',sep=""),
#     format='GTiff',overwrite=TRUE)
#   
#   lswi <- (band_vals[,2]/10000 - band_vals[,3]/10000)/(band_vals[,2]/10000 + band_vals[,3]/10000)
#   lswi_map <- setValues(nir,round(lswi*10000))
#   writeRaster(lswi_map,filename=paste(in_dirs[i],'/lswi',sep=""),
#     format='GTiff',overwrite=TRUE)
#   
#   sipi <- (band_vals[,3]/10000 - band_vals[,2]/10000)/(band_vals[,3]/10000 - band_vals[,1]/10000)
#   sipi_map <- setValues(nir,round(sipi*10000))
#   writeRaster(sipi_map,filename=paste(in_dirs[i],'/sipi',sep=""),
#     format='GTiff',overwrite=TRUE)
  
  mir <- (band_vals[,4]/10000)/(band_vals[,5]/10000)
  mir_map <- setValues(nir,round(mir*1000))
  writeRaster(mir_map,filename=paste(in_dirs[i],'/mir',sep=""),
    format='GTiff',overwrite=TRUE)
  
  count <- i
}
