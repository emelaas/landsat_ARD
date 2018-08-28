require(rgdal)
require(raster)
require(RColorBrewer)
require(classInt)
require(zyp)
require(Kendall)
library(scales)

can_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/Canada/','Canada')
can_shp_proj <- spTransform(can_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_shp <- readOGR('/projectnb/modislc/users/emelaas/scratch32/DAYMET/USA_Boundary/states_21basic/states.dbf','states')
usa_shp_proj <- spTransform(usa_shp, 
  CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
usa_crop <- crop(usa_shp_proj,extent(194214,2426712,-1660235,946090))
can_crop <- crop(can_shp_proj,extent(194214,2426712,-1660235,946090))


setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
RMSE4 <- raster('RMSE4.tif')
pval4 <- raster('pval4.tif')
cor4 <- raster('cor4.tif')
slope4 <- raster('slope_sig4.tif')

RMSE5 <- raster('RMSE5.tif')
pval5 <- raster('pval5.tif')
cor5 <- raster('cor5.tif')
slope5 <- raster('slope_sig5.tif')

RMSE6 <- raster('RMSE6.tif')
pval6 <- raster('pval6.tif')
cor6 <- raster('cor6.tif')
slope6 <- raster('slope_sig6.tif')

RMSE7 <- raster('RMSE7.tif')
pval7 <- raster('pval7.tif')
cor7 <- raster('cor7.tif')
slope7 <- raster('slope_sig7.tif')

nobs <- raster('sprLTM.tif')
nobs[nobs>0] <- 1

goodpix4 <- RMSE4
goodpix4[goodpix4>-9999] <- NA
goodpix4[RMSE4<5 & pval4<0.05 & cor4>0.5 & slope4==1] <- 1

goodpix5 <- RMSE5
goodpix5[goodpix5>-9999] <- NA
goodpix5[RMSE5<5 & pval5<0.05 & cor5>0.5 & slope5==1] <- 1

goodpix6 <- RMSE6
goodpix6[goodpix6>-9999] <- NA
goodpix6[RMSE6<5 & pval6<0.05 & cor6>0.5 & slope6==1] <- 1

goodpix7 <- RMSE7
goodpix7[goodpix7>-9999] <- NA
goodpix7[RMSE7<5 & pval7<0.05 & cor7>0.5 & slope7==1] <- 1

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/CONUS_ARD_grid/')
ard_shp <- readOGR(getwd(),'conus_ard_grid')
ard_shp_reproj <- spTransform(ard_shp,
  CRS('+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'))

tot_pix <- matrix(NA,length(ard_shp),5)
for (i in 1:length(ard_shp)){
  print(i)
  
  tile <- ard_shp_reproj[i,]
  tile_extr <- extract(stack(goodpix4,goodpix5,goodpix6,goodpix7,nobs),tile)

  if (sapply(tile_extr,is.null)==1){
    
  } else if (length(tile_extr[[1]][,1])==1){
    
  } else {
    tot_pix[i,] <- colSums(tile_extr[[1]],na.rm=TRUE)
  }
}

good_pix_frac <- matrix(NA,length(ard_shp),4)
for (i in 1:4){
  good_pix_frac[,i] <- tot_pix[,i]/tot_pix[,5]
}

w <- which(tot_pix[,5]>0)
good_pix_frac2 <- good_pix_frac[w,]
ard_shp2 <- ard_shp_reproj[w,]

#Generate map showing percentage of 500 m patches with significant spr/aut trends in each sidelap
cuts_s1 <- classIntervals(good_pix_frac2[,1],dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
plotclr <- rev(brewer.pal(9,'PuOr'))
colcode_s1 <- findColours(cuts_s1,plotclr)

cuts_s2 <- classIntervals(good_pix_frac2[,2],dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
plotclr <- rev(brewer.pal(9,'PuOr'))
colcode_s2 <- findColours(cuts_s2,plotclr)

cuts_s3 <- classIntervals(good_pix_frac2[,3],dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
plotclr <- rev(brewer.pal(9,'PuOr'))
colcode_s3 <- findColours(cuts_s3,plotclr)

cuts_s4 <- classIntervals(good_pix_frac2[,4],dataPrecision=3,style='fixed',
  fixedBreaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
plotclr <- rev(brewer.pal(9,'PuOr'))
colcode_s4 <- findColours(cuts_s4,plotclr)

#x11(h=8,w=8)
pdf(h=8,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/ard_tile_goodpix_opt.pdf')
par(mfrow=c(2,2),mar=c(1,1,1,1))
plot(usa_crop,col='white',main='SW Jan 1 (Mean AGDD)')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray30',add=TRUE)
plot(can_crop,col='gray30',add=TRUE)
plot(ard_shp2,pch=21,cex=1.6,col=alpha(colcode_s1,0.8),bg=colcode_s1,add=TRUE,lwd=0.5)
mtext('a',side=3,line=0,cex=1,adj=0,font=2)
legend(1926712,-600000,title='',cex=0.685,legend=names(attr(colcode_s1,'table')),
  fill=attr(colcode_s1,'palette'),bty='n')

plot(usa_crop,col='white',main='SW Feb 1 (Mean AGDD)')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray30',add=TRUE)
plot(can_crop,col='gray30',add=TRUE)
plot(ard_shp2,pch=21,cex=1.6,col=alpha(colcode_s2,0.8),bg=colcode_s2,add=TRUE,lwd=0.5)
mtext('b',side=3,line=0,cex=1,adj=0,font=2)

plot(usa_crop,col='white',main='SW Mar 17 (Mean AGDD)')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray30',add=TRUE)
plot(can_crop,col='gray30',add=TRUE)
plot(ard_shp2,pch=21,cex=1.6,col=alpha(colcode_s3,0.8),bg=colcode_s3,add=TRUE,lwd=0.5)
mtext('c',side=3,line=0,cex=1,adj=0,font=2)

plot(usa_crop,col='white',main='Photo (Exp=6)')
rect(194214,-1660235,2426712,946090,col='lightblue')
plot(usa_crop,col='gray30',add=TRUE)
plot(can_crop,col='gray30',add=TRUE)
plot(ard_shp2,pch=21,cex=1.6,col=alpha(colcode_s4,0.8),bg=colcode_s4,add=TRUE,lwd=0.5,)
mtext('d',side=3,line=0,cex=1,adj=0,font=2)

dev.off()

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
below5 <- raster('below5_sum.tif')
nobs <- raster('nobs.tif')
below5_vals <- getValues(below5) 
nobs_vals <- getValues(nobs)
below5_vals[is.na(nobs_vals)==1] <- NA
below5_vals[below5_vals==7] <- 8

all_vals <- nobs_vals
all_vals[all_vals>0] <- 1
all_vals <- cbind(all_vals,all_vals,all_vals,all_vals)

s <- stack(goodpix4,goodpix5,goodpix6,goodpix7)
s_vals <- getValues(s)
 
goodpix.group <- matrix(NA,8,4)
tot.group <- matrix(NA,8,4)
for (i in 1:4){
  goodpix.group[,i] <- aggregate(s_vals[,i],by=list(below5_vals),FUN=sum,na.rm=TRUE)[,2]
  tot.group[,i] <- aggregate(all_vals[,i],by=list(below5_vals),FUN=sum,na.rm=TRUE)[,2]
}

perc.group <- goodpix.group/tot.group

models <- c('SW Jan 1 (Mean AGDD)','SW Feb 1 (Mean AGDD)','SW Mar 17 (Mean AGDD)','Photo (Exp=6)')
colors <- c('royalblue','springgreen','hotpink','orange')

#x11(h=5,w=8)
pdf(h=5,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/goodpix_opt_by_Tavg.pdf')
barCenters <- barplot(t(perc.group),beside=TRUE,names=c(0,1,2,3,4,5,6,7),ylim=c(0,1),
  main='% "Good" Grid Cells',
  col=colors, ylab='',xlab=expression('Months with T'[avg]*' < 5'~degree*'C'))
legend('topleft',legend=models,
  cex=0.8,fill=colors,
  ncol=1)
dev.off()
