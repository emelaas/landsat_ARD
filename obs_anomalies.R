library(ncdf4)
require(rgdal)
library(raster)

can_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/Canada/','Canada')
can_shp_proj <- spTransform(can_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
usa_shp_proj <- spTransform(usa_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_crop <- crop(usa_shp_proj,extent(194214,2426712,-1660235,946090))
can_crop <- crop(can_shp_proj,extent(194214,2426712,-1660235,946090))

my.colors2 = colorRampPalette(c("white","green"))

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
all_yrs <- as.character(matrix(NA,34,1))
for (year in 1984:2017){
  in_dirs_tile <- list.files(path=getwd(),
    pattern=glob2rx(paste("spr*",year,"*tif",sep="")),
    full.names=T,include.dirs=T,recursive=TRUE)
  all_yrs[(year-1983)] <- in_dirs_tile[1]
}

all_yrs_stack <- stack(all_yrs)
all_yrs_vals <- getValues(all_yrs_stack)

o_mean <- rowMeans(all_yrs_vals,na.rm=TRUE)
o_sd <- apply(all_yrs_vals,1,sd,na.rm=TRUE)
writeRaster(setValues(all_yrs_stack[[1]],o_sd),
  filename='/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/obsSD',
  format='GTiff',overwrite=TRUE)

o_mean_all <- replicate(34, o_mean)
o_sd_all <- replicate(34, o_sd)
o_Zscore <- (all_yrs_vals-o_mean_all)/o_sd_all