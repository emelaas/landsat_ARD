##This script will generate long-term mean and annual 32km maps of spring and autumn phenology
##for a given Landast overlap scene

system.time({
  library(ncdf4)
  require(rgdal)
  library(raster)
  library(foreach)
  library(iterators)
  library(doParallel)
  
  #Register the parallel backend
  registerDoParallel(16)
  
  args = commandArgs(trailingOnly=T)
  tile_name = args[1]
  #tile_name <- "h18v16"
  
  #Load observed, 1km data
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=''))
  load(file = "landsat2daymet")
  
  #Native Daymet projection
  daymet_crs <- CRS("+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +a=6378137 +rf=298.257223563 +lat_1=25 +lat_2=60")
  
  #Import overlap shapefile
  o <- readOGR(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/SHP/',sep=""),tile_name)
  
  #Reproject overlap shapefile
  o_reproj <- spTransform(o,daymet_crs)
  
  tmean_all <- foreach(yr = 1982:2017, .combine = rbind) %dopar% {
    
    print(yr)
    
    setwd('/projectnb/modislc/data/daymet')
    tmax <- brick(paste('daymet_v3_tmax_',yr,'_na.nc4',sep=''),var='tmax')
    tmin <- brick(paste('daymet_v3_tmin_',yr,'_na.nc4',sep=''),var='tmin')
    
    tmax.crop <- crop(tmax,extent(o_reproj))
    tmin.crop <- crop(tmin,extent(o_reproj))
    
    tmean <- getValues((tmax.crop+tmin.crop)/2)
    
    time <- seq(as.Date(paste(yr,"/1/1",sep='')), as.Date(paste(yr,"/12/31",sep='')), "days")
    time <- time[1:365]
    month <- as.numeric(substr(time,6,7))
    year <- as.numeric(substr(time,1,4))
    
    tmean <- data.frame(year,month,t(tmean))        
  }
  
  month <- tmean_all[,2]
  
  monthly_tmean <- aggregate(tmean_all,list(month),FUN=mean)
  monthly_tmean <- monthly_tmean[,-c(1:3)]
  monthly_tmean[monthly_tmean < 5] <- 1
  monthly_tmean[monthly_tmean >= 5] <- 0
  
  below5_sum <- as.numeric(apply(monthly_tmean,2,sum))
  
  group <- matrix(NA,length(below5_sum),1)
  group[below5_sum < 3] <- 1
  group[below5_sum >= 3 & below5_sum < 6] <- 2
  group[below5_sum >=6] <- 3
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=''))
  save(group,below5_sum,file = 'below5_sum')  
})
