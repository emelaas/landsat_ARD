##This script will generate long-term mean and annual 32km maps of spring and autumn phenology
##for a given Landast overlap scene

library(ncdf4)
require(rgdal)
library(raster)
library(foreach)
library(iterators)
library(lmodel2)
library(doParallel)

#Register the parallel backend
registerDoParallel(16)

args = commandArgs(trailingOnly=T)
tile_name = args[1]
#tile_name <- "h22v10"

# Loop through each model and each daymet tile
for (m in 6:6){
  
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
  
  pred.SPR <- foreach(yr = 1982:2017, .combine = rbind) %dopar% {
    
    print(yr)
    
    setwd('/projectnb/modislc/data/climate/daymet')
    tmax <- brick(paste('daymet_v3_tmax_',yr,'_na.nc4',sep=''),var='tmax')
    tmin <- brick(paste('daymet_v3_tmin_',yr,'_na.nc4',sep=''),var='tmin')
    dayl <- brick(paste('daymet_v3_dayl_',yr,'_na.nc4',sep=''),var='dayl')
    
    tmax.crop <- crop(tmax,extent(o_reproj))
    tmin.crop <- crop(tmin,extent(o_reproj))
    dayl.crop <- crop(dayl,extent(o_reproj))
    
    tmean <- getValues((tmax.crop+tmin.crop)/2)
    photo <- getValues(dayl.crop)/3660
    
    # (AT_Jan1_5C) Alternating model with Jan 1 start date, 5C base temp
    # Parameters: a, b, c
    if (m == 1){
      #par <- rbind(c(0.85, 760.32, -0.017))
      par <- rbind(c(6.00, 797.53, -0.015))
      
      R_f <- tmean - 5
      R_f[R_f < 0] <- 0
      S_f <- t(apply(R_f, 1, cumsum))
      
      R_c <- tmean - 5
      R_c[R_c >= 0] <- 0
      R_c[R_c < 0] <- 1
      S_c <- t(apply(R_c, 1, cumsum))
      
      diff <- abs(S_f - (par[1] + par[2] * exp(par[3] * S_c)))
      pSPR <- as.numeric(apply(diff, 1, which.min))
      
      # (M1s_Jan1_5C) Photoperiod model with January 1 start date, 5C base temp
      # Parameters: EXP, F*
    } else if (m == 2){
      #par <- rbind(c(7.69, 1512.74))
      par <- rbind(c(6.19, 1383.35))
      
      R_f <- tmean - 5
      R_f[R_f < 0] <- 0
      R_f <- R_f*(photo/10)^par[1]
      S_f <- t(apply(R_f,1,cumsum))
      diff <- abs(S_f - par[2])
      
      pSPR <- as.numeric(apply(diff, 1, which.min))
      
      #         # (TT_Jan1_5C) Thermal time model with January 1 start date, 5C base temp
      #         # Parameters: F*
      #       } else if (m == 3){
      #         #par <- rbind(c(212.21))
      #         par <- rbind(c(285.2155))
      # 
      #         R_f <- tmean - 5
      #         R_f[R_f < 0] <- 0
      #         S_f <- t(apply(R_f,1,cumsum))
      #         diff <- abs(S_f - par[1])
      # 
      #         pSPR <- as.numeric(apply(diff, 1, which.min))
      
      # (TT) Thermal time model with flexible start date/base temp
      # Parameters: t0, Tbase, F*
    } else if (m == 3){
      #par <- rbind(c(73, -3.17, 487.95))
      par <- rbind(c(76, -1.22, 472.58))
      
      R_f <- tmean - par[2]
      R_f[R_f < 0] <- 0
      R_f[,1:par[1]] <- 0
      S_f <- t(apply(R_f,1,cumsum))
      diff <- abs(S_f - par[3])
      
      pSPR <- as.numeric(apply(diff, 1, which.min))
      
      # Use mean AGDD from Jan 1 (0 deg base temperature)
    } else if (m == 4){
      
      setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=''))
      load(file = "daymet_AGDD_Jan1_SOS")
      o_mean_mat <- replicate(365,o_mean)
      
      R_f <- tmean - 0
      R_f[R_f < 0] <- 0
      S_f <- t(apply(R_f, 1, cumsum))
      diff <- abs(S_f - o_mean_mat)
      
      pSPR <- as.numeric(apply(diff, 1, which.min))
      
      # Use mean AGDD from Feb 1 (0 deg base temperature)
    } else if (m == 5){
      
      setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=''))
      load(file = "daymet_AGDD_Feb1_SOS")
      o_mean_mat <- replicate(365,o_mean)
      
      R_f <- tmean - 0
      R_f[R_f < 0] <- 0
      R_f[,1:31] <- 0
      S_f <- t(apply(R_f, 1, cumsum))
      diff <- abs(S_f - o_mean_mat)
      
      pSPR <- as.numeric(apply(diff, 1, which.min))  
      
      # Use mean AGDD from Mar 17 (0 deg base temperature)
    } else if (m == 6){
      
      setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=''))
      load(file = "daymet_AGDD_Mar17_SOS")
      o_mean_mat <- replicate(365,o_mean)
      
      R_f <- tmean - 0
      R_f[R_f < 0] <- 0
      R_f[,1:75] <- 0
      S_f <- t(apply(R_f, 1, cumsum))
      diff <- abs(S_f - o_mean_mat)
      
      pSPR <- as.numeric(apply(diff, 1, which.min))  
      
      # Use mean AGDD from Jan 1 with Photoperiod exp=6 (0 deg base temperature)
    } else if (m == 7){
      
      setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=''))
      load(file = "daymet_AGDD_Jan1photo_exp6_SOS")
      o_mean_mat <- replicate(365,o_mean)
      
      R_f <- tmean - 0
      R_f[R_f < 0] <- 0
      R_f <- R_f*(photo/10)^6
      S_f <- t(apply(R_f, 1, cumsum))
      diff <- abs(S_f - o_mean_mat)
      
      pSPR <- as.numeric(apply(diff, 1, which.min))        
    }
  }
  
  pred.SPR[pred.SPR==1] <- NA
  
  assign('all.pSPR',t(pred.SPR))
  assign('all.oSPR',obs.SPR[,4:39])
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
  save(all.oSPR,all.pSPR,file = paste("daymet_predict_v5_",m,sep=""))
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
  load(paste("daymet_predict_v5_",m,sep=""))
  load('daymet_AGDD')
  
  #Calculate statistics for anomalously early observations
  o_mean <- apply(all.AGDD,1,mean,na.rm=TRUE)
  o_sd <- apply(all.AGDD,1,sd,na.rm=TRUE)
  o_mean_all <- replicate(36, o_mean)
  o_sd_all <- replicate(36, o_sd)
  o_Zscore <- (all.AGDD-o_mean_all)/o_sd_all
  
  # calculate RMSE and correlation for each Daymet grid cell
  nobs.pix <- matrix(NA,nrow(all.oSPR),1)
  rmse.pix <- matrix(NA,nrow(all.oSPR),1)
  rmse_2sd.pix <- matrix(NA,nrow(all.oSPR),1)
  mbe.pix <- matrix(NA,nrow(all.oSPR),1)
  cor.pix <- matrix(NA,nrow(all.oSPR),2)
  slope.pix <- matrix(NA,nrow(all.oSPR),1)
  slope_sig.pix <- matrix(NA,nrow(all.oSPR),1)
  for (i in 1:nrow(all.oSPR)){
    print(c(m,i))
    
    nobs <- sum(is.na(all.oSPR[i,])==0 & is.na(all.pSPR[i,])==0)
    nobs.pix[i] <- nobs
    
    if (nobs > 5){
      rmse.pix[i] <- sqrt(mean((all.oSPR[i,]-all.pSPR[i,])^2,na.rm=TRUE))
      mbe.pix[i] <- mean(all.oSPR[i,]-all.pSPR[i,],na.rm=TRUE)
      
      corr <- cor(all.oSPR[i,],all.pSPR[i,],use='pairwise.complete.obs')
      corr.test <- cor.test(all.oSPR[i,],all.pSPR[i,])
      cor.pix[i,] <- c(corr,corr.test$p.value)
      
      y1 <- all.oSPR[i,]
      x1 <- all.pSPR[i,]
      w <- which(is.na(y1)==0 & is.na(x1)==0)
      
      if (length(w)>10){
        lm2 <- lmodel2(y1[w]~x1[w],nperm=1,range.y="relative",range.x="relative")
        lb <- lm2$confidence.intervals[3,4]
        ub <- lm2$confidence.intervals[3,5]
        if (is.na(lb)==0 & lb < 1 & ub > 1) slope_sig.pix[i] <- 1
        slope.pix[i] <- lm2$regression.results[3,3]
      }

      #slope.pix[i] <- sd(all.oSPR[i,],na.rm=TRUE)/sd(all.pSPR[i,],na.rm=TRUE)
      
      w <- which(o_Zscore[i,] > 2)
      rmse_2sd.pix[i] <- sqrt(mean((all.oSPR[i,w]-all.pSPR[i,w])^2,na.rm=TRUE))
    }
  }
   
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
  save(all.oSPR,all.pSPR,o_Zscore,
    rmse.pix,mbe.pix,cor.pix,slope.pix,slope_sig.pix,rmse_2sd.pix,nobs.pix,
    file = paste("daymet_predict_",m,sep=""))
}


