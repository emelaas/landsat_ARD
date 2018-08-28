library(ncdf4)
require(rgdal)
library(raster)
library(lmodel2)
library(plyr)

sidelaps <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/sidelap_mask.tif')

tiles <- read.table('~/Code/GitHub/ard_tiles.txt')
tiles <- as.character(as.matrix(droplevels(tiles)))

all.tiles.pSPR1 <- matrix(NA,1,36)
all.tiles.pSPR2 <- matrix(NA,1,36)
all.tiles.pSPR3 <- matrix(NA,1,36)
all.tiles.pSPR4 <- matrix(NA,1,36)
all.tiles.oSPR <- matrix(NA,1,36)
all.tiles.below5 <- matrix(NA,1,1)
for (i in 1:length(tiles)){
  print(i)
  
  tile_name <- tiles[i]
  
  #Import overlap shapefile
  o <- readOGR(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/SHP/',sep=""),tile_name)
  o_reproj <- spTransform(o,'+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
  sidelap.crop <- crop(sidelaps,o_reproj)
  sidelap.crop_vals <- getValues(sidelap.crop)
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
  load('daymet_predict_v5_1')
  all.tiles.oSPR <- rbind(all.tiles.oSPR,all.oSPR)
  all.tiles.pSPR1 <- rbind(all.tiles.pSPR1,all.pSPR)
  
  load('daymet_predict_v5_2')
  all.tiles.pSPR2 <- rbind(all.tiles.pSPR2,all.pSPR)
  
  load('daymet_predict_v5_3')
  all.tiles.pSPR3 <- rbind(all.tiles.pSPR3,all.pSPR)
  
  load('daymet_predict_v5_4')
  all.tiles.pSPR4 <- rbind(all.tiles.pSPR4,all.pSPR)
  
  load('below5_sum')
  below5_sum[sidelap.crop_vals==0] <- NA
  all.tiles.below5 <- rbind(all.tiles.below5,as.matrix(below5_sum))
}

w1 <- which(all.tiles.below5<3)
x11 <- all.tiles.pSPR1[w1,]
y1 <- all.tiles.oSPR[w1,]
x21 <- all.tiles.pSPR2[w1,]
x31 <- all.tiles.pSPR3[w1,]
x41 <- all.tiles.pSPR4[w1,]

w2 <- which(all.tiles.below5>=3 & all.tiles.below5<=5)
x12 <- all.tiles.pSPR1[w2,]
y2 <- all.tiles.oSPR[w2,]
x22 <- all.tiles.pSPR2[w2,]
x32 <- all.tiles.pSPR3[w2,]
x42 <- all.tiles.pSPR4[w2,]

w3 <- which(all.tiles.below5>=6)
x13 <- all.tiles.pSPR1[w3,]
y3 <- all.tiles.oSPR[w3,]
x23 <- all.tiles.pSPR2[w3,]
x33 <- all.tiles.pSPR3[w3,]
x43 <- all.tiles.pSPR4[w3,]

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/R_files/')
save(x11,x21,x31,x41,x12,x22,x32,x42,x13,x23,x33,x43,y1,y2,y3,file='scatter_mean_below5')

y1mean <- colMeans(y1,na.rm=TRUE)
x11mean <- colMeans(x11,na.rm=TRUE)
x21mean <- colMeans(x21,na.rm=TRUE)
x31mean <- colMeans(x31,na.rm=TRUE)
x41mean <- colMeans(x41,na.rm=TRUE)

y2mean <- colMeans(y2,na.rm=TRUE)
x12mean <- colMeans(x12,na.rm=TRUE)
x22mean <- colMeans(x22,na.rm=TRUE)
x32mean <- colMeans(x32,na.rm=TRUE)
x42mean <- colMeans(x42,na.rm=TRUE)

y3mean <- colMeans(y3,na.rm=TRUE)
x13mean <- colMeans(x13,na.rm=TRUE)
x23mean <- colMeans(x23,na.rm=TRUE)
x33mean <- colMeans(x33,na.rm=TRUE)
x43mean <- colMeans(x43,na.rm=TRUE)

#x11(h=9,w=12)
pdf(h=9,w=12,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/scatter_mean_below5.pdf')
par(mfrow=c(3,4),mar=c(4,4,1.5,1))

plot(x13mean,y3mean,pch=21,cex=1.5,col='black',bg='royalblue',xlim=c(130,160),ylim=c(130,160),
  ylab='Observed SOS (DOY)',xlab='Predicted SOS (DOY)',main='Chill')
mtext('a',side=3,line=0,cex=1,adj=-0.2,font=2)
text(138,156,expression('6-8 \nMonths T'[avg]*' < 5'~degree*'C'),cex=1.1)
lm2 <- lmodel2(y3mean~x13mean)
rmse <- sqrt(mean((y3mean-x13mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x23mean,y3mean,pch=21,cex=1.5,col='black',bg='springgreen',xlim=c(130,160),ylim=c(130,160),
  ylab='',xlab='',main='Photo')
mtext('b',side=3,line=0,cex=1,adj=-0.2,font=2)
lm2 <- lmodel2(y3mean~x23mean)
rmse <- sqrt(mean((y3mean-x23mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x33mean,y3mean,pch=21,cex=1.5,col='black',bg='hotpink',xlim=c(130,160),ylim=c(130,160),
  ylab='',xlab='',main='SW Mar 17')
mtext('c',side=3,line=0,cex=1,adj=-0.2,font=2)
lm2 <- lmodel2(y3mean~x33mean)
rmse <- sqrt(mean((y3mean-x33mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x43mean,y3mean,pch=21,cex=1.5,col='black',bg='orange',xlim=c(130,160),ylim=c(130,160),
  ylab='',xlab='',main='SW Jan 1 (mean AGDD)')
mtext('d',side=3,line=0,cex=1,adj=-0.2,font=2)
lm2 <- lmodel2(y3mean~x43mean)
rmse <- sqrt(mean((y3mean-x43mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x12mean,y2mean,pch=21,cex=1.5,col='black',bg='royalblue',xlim=c(105,145),ylim=c(105,145),
  ylab='',xlab='')
mtext('e',side=3,line=0,cex=1,adj=-0.2,font=2)
text(115,141,expression('3-5 \nMonths T'[avg]*' < 5'~degree*'C'),cex=1.1)
lm2 <- lmodel2(y2mean~x12mean)
rmse <- sqrt(mean((y2mean-x12mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x22mean,y2mean,pch=21,cex=1.5,col='black',bg='springgreen',xlim=c(105,145),ylim=c(105,145),
  ylab='',xlab='')
mtext('f',side=3,line=0,cex=1,adj=-0.2,font=2)
lm2 <- lmodel2(y2mean~x22mean)
rmse <- sqrt(mean((y2mean-x22mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x32mean,y2mean,pch=21,cex=1.5,col='black',bg='hotpink',xlim=c(105,145),ylim=c(105,145),
  ylab='',xlab='')
mtext('g',side=3,line=0,cex=1,adj=-0.2,font=2)
lm2 <- lmodel2(y2mean~x32mean)
rmse <- sqrt(mean((y2mean-x32mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x42mean,y2mean,pch=21,cex=1.5,col='black',bg='orange',xlim=c(105,145),ylim=c(105,145),
  ylab='',xlab='')
mtext('h',side=3,line=0,cex=1,adj=-0.2,font=2)
lm2 <- lmodel2(y2mean~x42mean)
rmse <- sqrt(mean((y2mean-x42mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x11mean,y1mean,pch=21,cex=1.5,col='black',bg='royalblue',xlim=c(85,120),ylim=c(85,120),
  ylab='',xlab='')
mtext('i',side=3,line=0,cex=1,adj=-0.2,font=2)
text(94,116,expression('0-2 \nMonths T'[avg]*' < 5'~degree*'C'),cex=1.1)
lm2 <- lmodel2(y1mean~x11mean)
rmse <- sqrt(mean((y1mean-x11mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x21mean,y1mean,pch=21,cex=1.5,col='black',bg='springgreen',xlim=c(85,120),ylim=c(85,120),
  ylab='',xlab='')
mtext('j',side=3,line=0,cex=1,adj=-0.2,font=2)
lm2 <- lmodel2(y1mean~x21mean)
rmse <- sqrt(mean((y1mean-x21mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x31mean,y1mean,pch=21,cex=1.5,col='black',bg='hotpink',xlim=c(85,120),ylim=c(85,120),
  ylab='',xlab='')
mtext('k',side=3,line=0,cex=1,adj=-0.2,font=2)
lm2 <- lmodel2(y1mean~x31mean)
rmse <- sqrt(mean((y1mean-x31mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

plot(x41mean,y1mean,pch=21,cex=1.5,col='black',bg='orange',xlim=c(85,120),ylim=c(85,120),
  ylab='',xlab='')
mtext('l',side=3,line=0,cex=1,adj=-0.2,font=2)
lm2 <- lmodel2(y1mean~x41mean)
rmse <- sqrt(mean((y1mean-x41mean)^2,na.rm=TRUE))
abline(lm2$regression.results[3,2],lm2$regression.results[3,3],lwd=0.75)
rp = vector('expression',4)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(lm2$rsquare,dig=2),MYVALUE2 = format(lm2$P.param,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse,dig=2)))[2]
rp[3] = substitute(expression(paste(beta[1] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,3], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,4], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,5], digits=2)))[2]
rp[4] = substitute(expression(paste(beta[0] == MYVALUE1,' (',MYVALUE2,',',MYVALUE3,')')),
  list(MYVALUE1 = format(lm2$regression.results[3,2], digits=2),
    MYVALUE2 = format(lm2$confidence.intervals[3,2], digits=2),
    MYVALUE3 = format(lm2$confidence.intervals[3,3], digits=2)))[2]
legend('bottomright', legend = rp, bty = 'n',cex=1.1)
abline(0,1,lty=2,lwd=0.5)

dev.off()