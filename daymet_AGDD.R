##This script will generate long-term mean and annual 32km maps of spring and autumn phenology
##for a given Landast overlap scene

library(ncdf4)
require(rgdal)
library(raster)

args = commandArgs(trailingOnly=T)
tile_name = args[1]

#tile_name <- "h18v16"

# Load in Daymet tile #s
setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/CLIM/',sep=''))
in_dirs <- list.files(path=getwd(),pattern=glob2rx("tmin_1982*"),full.names=T,include.dirs=T)
tiles <- as.numeric(substr(in_dirs,nchar(in_dirs[1])-7,nchar(in_dirs[1])-3))

# Loop through each model and each daymet tile

for (i in 1:length(tiles)){

  #Load observed, 1km data
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
  load(file = paste("landsat2daymet_",tiles[i],sep=''))

  lon <- raster(in_dirs[i],varname='lon')
  lat <- raster(in_dirs[i],varname='lat')
  cells <- setValues(lon,seq(1,ncell(lon)))

  #Generate lat/lon points from rasters and reproject to native Daymet projection
  plat <- rasterToPoints(lat)
  plon <- rasterToPoints(lon)

  lonlat <- cbind(plon[,3], plat[,3])
  lonlat <- SpatialPoints(lonlat, proj4string = CRS("+proj=longlat +datum=WGS84"))

  #Native Daymet projection
  daymet_crs <- CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  plonlat <- spTransform(lonlat, CRSobj = daymet_crs)

  projection(lat) <- daymet_crs
  extent(lat) <- extent(plonlat)
  projection(lon) <- daymet_crs
  extent(lon) <- extent(plonlat)
  projection(cells) <- daymet_crs
  extent(cells) <- extent(plonlat)

  #Import overlap shapefile
  o <- readOGR(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/SHP/',sep=""),tile_name)

  #Reproject overlap shapefile
  o_reproj <- spTransform(o,daymet_crs)

  lat.crop <- try(crop(lat,extent(o_reproj)))
  if(!inherits(lat.crop,"try-error")){
    lat_vals = getValues(lat.crop)

    lon.crop <- crop(lon,extent(o_reproj))
    lon_vals = getValues(lon.crop)

    cells.crop <- crop(cells,extent(o_reproj))
    cells_vals = getValues(cells.crop)

    AGDD <- matrix(NA,nrow(obs.SPR),36)
    for (yr in 1982:2017){

      print(yr)

      setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/CLIM/',sep=''))
      tmax <- brick(paste('tmax_',yr,'_',tiles[i],'.nc',sep=''),var='tmax')
      tmin <- brick(paste('tmin_',yr,'_',tiles[i],'.nc',sep=''),var='tmin')
      dayl <- brick(paste('dayl_',yr,'_',tiles[i],'.nc',sep=''),var='dayl')

      projection(tmax) <- daymet_crs
      extent(tmax) <- extent(plonlat)
      projection(tmin) <- daymet_crs
      extent(tmin) <- extent(plonlat)
      projection(dayl) <- daymet_crs
      extent(dayl) <- extent(plonlat)

      tmax.crop <- crop(tmax,extent(o_reproj))
      tmin.crop <- crop(tmin,extent(o_reproj))
      dayl.crop <- crop(dayl,extent(o_reproj))

      tmean <- getValues((tmax.crop+tmin.crop)/2)
      tmean[tmean<0] <- 0

      for (j in 1:nrow(obs.SPR)){
        if (is.na(obs.SPR[j,2])==0 & obs.SPR[j,2]>60){
          AGDD[j,(yr-1981)] <- sum(tmean[j,(obs.SPR[j,2]-60):obs.SPR[j,2]])
        }
      }
    }

  } else {
    AGDD <- matrix(NA,1,36)
  }

  if (i == 1){
    assign('all.AGDD',AGDD)
  } else {
    all.AGDD <- rbind(all.AGDD,AGDD)
  }
}

setwd(paste('/projectnb/modislc/projects/landsat_senintel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
save(all.AGDD,file = "daymet_AGDD")
