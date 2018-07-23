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

all_yrs_vals[is.na(all_yrs_vals) == 1] <- 9999
all_yrs_vals <- -1*all_yrs_vals
early_yr <- max.col(all_yrs_vals)

count <- all_yrs_vals
count[count!=(-9999)] <- 0
count[count==(-9999)] <- 1
count_sum <- rowSums(count)

early_yr[count_sum>=30] <- NA
early_yr_map <- setValues(all_yrs_stack[[1]],early_yr)

for (i in 1984:2017){
  jpeg(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/early_yr/early_',i,'.jpeg',sep=''))
  plot(usa_crop,col='white',main=i)
  rect(194214,-1660235,2426712,946090,col='lightblue')
  plot(usa_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(can_crop,col='gray50',add=TRUE,lwd=0.5)
  plot(early_yr_map==(i-1983),legend=FALSE,add=TRUE,col=my.colors2(2))
  dev.off()
}
