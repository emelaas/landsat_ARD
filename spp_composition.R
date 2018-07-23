library(raster)
library(rgdal)
library(gdalUtils)
library(foreach)
library(iterators)
library(doParallel)
library(ncdf4)

args = commandArgs(trailingOnly=T)
tile = args[1]
#tile <- 'h18v06'

print(tile)

#Register the parallel backend
registerDoParallel(8)

# Load lat/lon coordinates for CLM
setwd('/projectnb/modislc/users/emelaas/MSB/CLM/')
ncin <- nc_open('cesm_LatLon_file_fv_0p45_0p625.nc')

# Clip extreme values (e.g., -180/0 lon, -90/90 lat)
lon <- ncvar_get(ncin,"lon")
lat <- ncvar_get(ncin,"lat")
lon <- lon[-c(1,289)]
lat <- lat[-c(1,383)]

# Generate matrix 
lon_mat <- matrix(lon,nrow=381,ncol=length(lon),byrow=TRUE)
lat_mat <- t(matrix(lat,nrow=574,ncol=length(lat),byrow=TRUE))
dim(lon_mat) <- c(381*574,1)
dim(lat_mat) <- c(381*574,1)

# Generate spatial polygons for aggregation and reproject to sinusoidal
t <- data.frame(lat_mat,lon_mat)
names(t) <- c('lat','lon')
coordinates(t) <- ~lon+lat
proj4string(t) <- CRS("+proj=longlat +datum=WGS84")
gridded(t) = TRUE
r = raster(t)
temp <- seq(1,ncell(r),1)
r <- setValues(r,temp)
r <- crop(r,extent(r,1,190,1,594))
r_poly <- rasterToPolygons(r)
r_reproj <- spTransform(r_poly,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))

setwd('/projectnb/modislc/users/emelaas/MSB/LAI/tile_means/halfdegree/')
in_dirs_complete <- list.files(path=getwd(),
  pattern=glob2rx(tile),full.names=T,include.dirs=T,recursive=T)

setwd(paste('/projectnb/modislc/users/emelaas/MSB/LAI/qa_mask/',tile,sep=''))
in_dirs_mod15 <- list.files(path=getwd(),
  pattern=glob2rx('lai*'),full.names=T,include.dirs=T,recursive=T)

complete_dates <- substring(in_dirs_complete,64,70)
all_dates <- substring(in_dirs_mod15,68,74)
todo_dates <- which(!all_dates %in% complete_dates)

if (length(todo_dates)>0) {
  
  print(todo_dates)
  
  all_tile_means <- foreach(i = todo_dates, .combine = rbind) %dopar% {
    
    #Get list of individual bands
    LAI <- raster(in_dirs_mod15[i])
    LAI <- LAI/10
    
    date <- substr(in_dirs_mod15[i],68,74)
    year <- substr(in_dirs_mod15[i],68,71)
    
    year[year<2001] <- 2001
    year[year>2016] <- 2016
    
    setwd('/projectnb/modislc/users/emelaas/MSB/MCD12Q1/')
    hdfs_lc <- list.files(path=getwd(),pattern=glob2rx(paste('*A',year,'*',tile,"*.hdf",sep='')),
      full.names=TRUE)
    sds <- get_subdatasets(hdfs_lc)
    LC <- raster(sds[1])
    proj4string(LC) <- CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs')
    LC[LC==0] <- NA
    LC[LC>20] <- NA
    
    lc_means <- matrix(NA,length(r_reproj),18)
    lc_frac <- matrix(NA,length(r_reproj),18)
    num_pix <- matrix(NA,length(r_reproj),1)
    for (j in 1:length(r_reproj)){
      #print(j)
      if (j%%10000==0) print(c(i,j))
      
      r_extr <- extract(stack(LAI,LC),r_reproj[j,])
      if (sapply(r_extr,is.null)==1){
        tile_mean <- NA
        
      } else if (length(r_extr[[1]][,1])==1){
        tile_mean <- NA
        
      } else {
        num_pix[j] <- length(which(is.na(r_extr[[1]][,1])==0))
        
        df <- as.data.frame(cbind(r_extr[[1]][,1],r_extr[[1]][,2]))
        colnames(df) <- c('lai','lc')
        tile_mean <- aggregate(df$lai,by=list(df$lc),FUN=mean,na.rm=TRUE)
        lc_means[j,tile_mean$Group.1] <- tile_mean$x
        lc_means[j,18] <- mean(df$lai,na.rm=TRUE)
        
        tile_lc_tot <- table(df$lc)
        tmp <- matrix(NA,1,18)
        lc_frac[j,as.numeric(names(tile_lc_tot))] <- tile_lc_tot/sum(tile_lc_tot)        
      }
    }
    
    # Save model calibration results
    save(lc_means,lc_frac,num_pix,
      file=paste('/projectnb/modislc/users/emelaas/MSB/LAI/tile_means/halfdegree/',
        date,'/',tile,sep = ''))
    
    tile_means <- i
    
  }  
}

