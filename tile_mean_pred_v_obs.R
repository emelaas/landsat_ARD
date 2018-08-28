tiles <- read.table('~/Code/GitHub/ard_tiles.txt')
tiles <- as.character(as.matrix(droplevels(tiles)))

for (i in 1:length(tiles)){
  print(i)
  
  tile_name <- tiles[i]
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
  pdf(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/tile_mean_pred_v_obs/',
    tile_name,'.pdf',sep=''))
  par(mfrow=c(3,3),mar=c(2,2,2,2))
  
  for(m in 1:7){
    load(paste('daymet_predict_v5_',m,sep=''))
    if(m==1){
      x <- colMeans(all.pSPR,na.rm=TRUE)
      y <- colMeans(all.oSPR,na.rm=TRUE)
      xlow <- min(c(x,y),na.rm=TRUE)-5
      xhi <- max(c(x,y),na.rm=TRUE)+5
      ylow <- min(c(x,y),na.rm=TRUE)-5
      yhi <- max(c(x,y),na.rm=TRUE)+5
      plot(x,y,pch=16,xlim=c(xlow,xhi),ylim=c(ylow,yhi))
      abline(0,1)
    } else {
      x1 <- colMeans(all.pSPR,na.rm=TRUE)
      y1 <- colMeans(all.oSPR,na.rm=TRUE)
      plot(x1,y1,pch=16,xlim=c(xlow,xhi),ylim=c(ylow,yhi))
      abline(0,1)
    }


  }
  dev.off()
}






tiles <- read.table('~/Code/GitHub/ard_tiles.txt')
tiles <- as.character(as.matrix(droplevels(tiles)))
nums <- na.omit(as.numeric(unlist(strsplit(unlist(tiles),"[^0-9]+"))))
dim(nums) <- c(2,length(tiles))
tiles <- tiles[which(nums[1,]==20)]

x11(h=length(tiles),w=5)
# pdf(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/figures/tile_mean_pred_v_obs/all_models_h',
#   20,'.pdf',sep=''))
par(mfrow=c(length(tiles),4),mar=c(0,0,0,0))
for (i in 1:length(tiles)){
  print(i)
  
  tile_name <- tiles[i]
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))

  
  for(m in 1:4){
    load(paste('daymet_predict_v5_',m,sep=''))
    if(m==1){
      x <- colMeans(all.pSPR,na.rm=TRUE)
      y <- colMeans(all.oSPR,na.rm=TRUE)
      xlow <- min(c(x,y),na.rm=TRUE)-5
      xhi <- max(c(x,y),na.rm=TRUE)+5
      ylow <- min(c(x,y),na.rm=TRUE)-5
      yhi <- max(c(x,y),na.rm=TRUE)+5
      plot(x,y,pch=16,xlim=c(xlow,xhi),ylim=c(ylow,yhi),xaxt='n',yaxt='n',xlab='n',ylab='n')
      abline(0,1,lty=2,lwd=0.5)
    } else {
      x1 <- colMeans(all.pSPR,na.rm=TRUE)
      y1 <- colMeans(all.oSPR,na.rm=TRUE)
      plot(x1,y1,pch=16,xlim=c(xlow,xhi),ylim=c(ylow,yhi),xaxt='n',yaxt='n',xlab='n',ylab='n')
      abline(0,1,lty=2,lwd=0.5)
    }
    
    
  }
}

dev.off()