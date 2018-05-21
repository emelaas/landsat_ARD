require(rgdal)
library(raster)

wr2 <- readOGR('/projectnb/modislc/projects/te_phenology/landsat_scenes/wrs2_descending','wrs2_descending')
wr2_proj <- spTransform(wr2, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("nobs*tif"),full.names=T,include.dirs=T,recursive=TRUE)
nobs <- raster(in_dirs_tile[1])
count_map <- setValues(nobs,seq(1,ncell(nobs)))
count <- matrix(0,ncell(nobs),1)

for (path in 30:10){
  for (row in 26:40){
    print(c(path,row))
    
    PR1 <- path*1000+row
    PR2 <- (path-1)*1000+row
    PR3 <- (path-1)*1000+(row+1)
    wr2_1 <- wr2_proj[which(wr2_proj$PR==PR1),]
    wr2_2 <- wr2_proj[which(wr2_proj$PR==PR2),]
    wr2_3 <- wr2_proj[which(wr2_proj$PR==PR3),]
    
    int1 <- intersect(wr2_1,wr2_2)
    int2 <- intersect(wr2_1,wr2_3)
    
    poly1 <- extract(count_map,int1)
    pixels <- unlist(poly1)
    count[pixels] <- count[pixels]+1
    
    poly2 <- extract(count_map,int2)
    pixels <- unlist(poly2)
    count[pixels] <- count[pixels]+1
  }
}

count <- setValues(count_map,count)
writeRaster(count,
  filename='/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sidelap_mask',
  format='GTiff',overwrite=TRUE)