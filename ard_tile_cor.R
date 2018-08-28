require(rgdal)
require(raster)
require(RColorBrewer)
require(classInt)
require(zyp)
require(Kendall)
library(scales)
library(foreach)
library(iterators)
library(doParallel)

#Register the parallel backend
registerDoParallel(16)

can_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/Canada/','Canada')
can_shp_proj <- spTransform(can_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
usa_shp_proj <- spTransform(usa_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_crop <- crop(usa_shp_proj,extent(194214,2426712,-1660235,946090))
can_crop <- crop(can_shp_proj,extent(194214,2426712,-1660235,946090))


setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
cor1 <- raster('cor1.tif')
cor2 <- raster('cor2.tif')
cor3 <- raster('cor3.tif')
cor4 <- raster('cor4.tif')

nobs <- raster('sprLTM.tif')
nobs[nobs>0] <- 1

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/CONUS_ARD_grid/')
ard_shp <- readOGR(getwd(),'conus_ard_grid')
ard_shp_reproj <- spTransform(ard_shp,
  CRS('+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'))

tot_pix_all <- foreach(i = 1:length(ard_shp), .combine = rbind) %dopar% {
  print(i)
  
  tile <- ard_shp_reproj[i,]
  
  tile_extr <- extract(nobs,tile)
  if (length(which(tile_extr[[1]]>0)>0)){
    tile_extr <- extract(stack(cor1,cor2,cor3,cor4,nobs),tile)
    
    if (sapply(tile_extr,is.null)==1){
      tot_pix <- as.numeric(matrix(NA,1,5))
    } else if (length(tile_extr[[1]][,1])==1){
      tot_pix <- as.numeric(matrix(NA,1,5))
    } else {
      tot_pix <- colMeans(tile_extr[[1]],na.rm=TRUE)
    }
  } else {
    tot_pix <- as.numeric(matrix(NA,1,5))
  }
}

save(tot_pix_all,file='/projectnb/modislc/projects/landsat_sentinel/ARD/R_files/tot_pix_all_cor')

w <- which(tot_pix_all[,5]>0)
tot_pix2 <- tot_pix_all[w,]
ard_shp2 <- ard_shp_reproj[w,]

#Generate map showing percentage of 500 m patches with significant spr/aut trends in each sidelap
cuts_s1 <- classIntervals(tot_pix2[,1],dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.5,0.6,0.7,0.8,0.9,1))
plotclr <- rev(brewer.pal(9,'YlOrRd'))
colcode_s1 <- findColours(cuts_s1,plotclr)

cuts_s2 <- classIntervals(tot_pix2[,2],dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.5,0.6,0.7,0.8,0.9,1))
plotclr <- rev(brewer.pal(9,'YlOrRd'))
colcode_s2 <- findColours(cuts_s2,plotclr)

cuts_s3 <- classIntervals(tot_pix2[,3],dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.5,0.6,0.7,0.8,0.9,1))
plotclr <- rev(brewer.pal(9,'YlOrRd'))
colcode_s3 <- findColours(cuts_s3,plotclr)

cuts_s4 <- classIntervals(tot_pix2[,4],dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.5,0.6,0.7,0.8,0.9,1))
plotclr <- rev(brewer.pal(9,'YlOrRd'))
colcode_s4 <- findColours(cuts_s4,plotclr)

#x11(h=8,w=8)
pdf(h=8,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/ard_tile_cor.pdf')
par(mfrow=c(2,2),mar=c(1,1,1,1))
plot(usa_crop,col='white',main='Chill')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray30',add=TRUE)
plot(can_crop,col='gray30',add=TRUE)
plot(ard_shp2,pch=21,cex=1.6,col=alpha(colcode_s1,0.8),bg=colcode_s1,add=TRUE,lwd=0.5)
mtext('a',side=3,line=0,cex=1,adj=0,font=2)
legend(1926712,-600000,title='',cex=0.685,legend=names(attr(colcode_s1,'table')),
  fill=attr(colcode_s1,'palette'),bty='n')

plot(usa_crop,col='white',main='Photo')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray30',add=TRUE)
plot(can_crop,col='gray30',add=TRUE)
plot(ard_shp2,pch=21,cex=1.6,col=alpha(colcode_s2,0.8),bg=colcode_s2,add=TRUE,lwd=0.5)
mtext('b',side=3,line=0,cex=1,adj=0,font=2)

plot(usa_crop,col='white',main='SW Mar 17')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray30',add=TRUE)
plot(can_crop,col='gray30',add=TRUE)
plot(ard_shp2,pch=21,cex=1.6,col=alpha(colcode_s3,0.8),bg=colcode_s3,add=TRUE,lwd=0.5)
mtext('c',side=3,line=0,cex=1,adj=0,font=2)

plot(usa_crop,col='white',main='SW Jan 1 (Mean AGDD)')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray30',add=TRUE)
plot(can_crop,col='gray30',add=TRUE)
plot(ard_shp2,pch=21,cex=1.6,col=alpha(colcode_s4,0.8),bg=colcode_s4,add=TRUE,lwd=0.5,)
mtext('d',side=3,line=0,cex=1,adj=0,font=2)

dev.off()
