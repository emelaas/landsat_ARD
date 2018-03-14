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

  # Load in Daymet tile #s
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/CLIM/',sep=''))
  in_dirs <- list.files(path=getwd(),pattern=glob2rx("tmin_1982*"),full.names=T,include.dirs=T)
  tiles <- as.numeric(substr(in_dirs,nchar(in_dirs[1])-7,nchar(in_dirs[1])-3))

  # Number of models tested
  num_mod <- 4

  # Initialize matrices for model performance stats
  mbe <- matrix(NA,num_mod,1)
  rmse <- matrix(NA,num_mod,1)
  cor <- matrix(NA,num_mod,1)
  mbe_2sd <- matrix(NA,num_mod,1)
  rmse_2sd <- matrix(NA,num_mod,1)
  cor_2sd <- matrix(NA,num_mod,1)
  good.pix <- matrix(NA,num_mod,1)

  # Loop through each model and each daymet tile
  for (m in 1:num_mod){
    for (i in 1:length(tiles)){

      #Load observed, 1km data
      setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=''))
      load(file = paste("landsat2daymet_",tiles[i],sep=''))

      lon <- raster(in_dirs[i],varname='lon')
      lat <- raster(in_dirs[i],varname='lat')
      cells <- setValues(lon,seq(1,ncell(lon)))

      #Generate lat/lon points from rasters and reproject to native Daymet projection
      plat <- rasterToPoints(lat)
      plon <- rasterToPoints(lon)

      lonlat <- cbind(plon[,3], plat[,3])
      lonlat <- SpatialPoints(lonlat, proj4string = CRS("+proj=longlat +datum=WGS84"))

      #Native Daymet projection
      daymet_crs <- CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
      plonlat <- spTransform(lonlat, CRSobj = daymet_crs)

      projection(lat) <- daymet_crs
      extent(lat) <- extent(plonlat)
      projection(lon) <- daymet_crs
      extent(lon) <- extent(plonlat)
      projection(cells) <- daymet_crs
      extent(cells) <- extent(plonlat)

      #Import overlap shapefile
      o <- readOGR(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
        tile_name,'/SHP/',sep=""),tile_name)

      #Reproject overlap shapefile
      o_reproj <- spTransform(o,daymet_crs)


      lat.crop <- try(crop(lat,extent(o_reproj)))
      if(!inherits(lat.crop,"try-error")){
        lat_vals = getValues(lat.crop)

        lon.crop <- crop(lon,extent(o_reproj))
        lon_vals = getValues(lon.crop)

        cells.crop <- crop(cells,extent(o_reproj))
        cells_vals = getValues(cells.crop)

        pred.SPR <- foreach(yr = 1982:2017, .combine = rbind) %dopar% {

          print(yr)

          setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/CLIM/',sep=''))
          tmax <- brick(paste('tmax_',yr,'_',tiles[i],'.nc',sep=''),var='tmax')
          tmin <- brick(paste('tmin_',yr,'_',tiles[i],'.nc',sep=''),var='tmin')
          dayl <- brick(paste('dayl_',yr,'_',tiles[i],'.nc',sep=''),var='dayl')

          projection(tmax) <- daymet_crs
          extent(tmax) <- extent(plonlat)
          projection(tmin) <- daymet_crs
          extent(tmin) <- extent(plonlat)
          projection(dayl) <- daymet_crs
          extent(dayl) <- extent(plonlat)

          tmax.crop <- crop(tmax,extent(o_reproj))
          tmin.crop <- crop(tmin,extent(o_reproj))
          dayl.crop <- crop(dayl,extent(o_reproj))

          tmean <- getValues((tmax.crop+tmin.crop)/2)
          photo <- getValues(dayl.crop)/3660

          # (AT_Jan1_5C) Alternating model with Jan 1 start date, 5C base temp
          # Parameters: a, b, c
          if (m == 1){
            par <- rbind(c(0.85, 760.32, -0.017))

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
            par <- rbind(c(7.69, 1512.74))

            R_f <- tmean - 5
            R_f[R_f < 0] <- 0
            R_f <- R_f*(photo/10)^par[1]
            S_f <- t(apply(R_f,1,cumsum))
            diff <- abs(S_f - par[2])

            pSPR <- as.numeric(apply(diff, 1, which.min))

            # (TT_Jan1_5C) Thermal time model with January 1 start date, 5C base temp
            # Parameters: F*
          } else if (m == 3){
            par <- rbind(c(212.21))

            R_f <- tmean - 5
            R_f[R_f < 0] <- 0
            S_f <- t(apply(R_f,1,cumsum))
            diff <- abs(S_f - par[1])

            pSPR <- as.numeric(apply(diff, 1, which.min))

            # (TT) Thermal time model with flexible start date/base temp
            # Parameters: t0, Tbase, F*
          } else if (m == 4){
            par <- rbind(c(73, -3.17, 487.95))

            R_f <- tmean - par[2]
            R_f[R_f < 0] <- 0
            R_f[,1:par[1]] <- 0
            S_f <- t(apply(R_f,1,cumsum))
            diff <- abs(S_f - par[3])

            pSPR <- as.numeric(apply(diff, 1, which.min))
          }
        }
      } else {
        pred.SPR <- matrix(NA,32,1)
      }

      pred.SPR[pred.SPR==1] <- NA

      if (i == 1){
        assign('all.pSPR',t(pred.SPR))
        assign('all.oSPR',obs.SPR[,4:39])
      } else {
        all.pSPR <- rbind(all.pSPR,t(pred.SPR))
        all.oSPR <- rbind(all.oSPR,obs.SPR[,4:39])
      }
    }

    setwd(paste('/projectnb/modislc/projects/landsat_senintel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
    save(all.oSPR,all.pSPR,file = paste("daymet_predict_",m,sep=""))

    load(paste("daymet_predict_",m,sep=""))
    load('daymet_AGDD')

    # calculate RMSE and correlation for each Daymet grid cell
    rmse.pix <- matrix(NA,nrow(all.oSPR),1)
    cor.pix <- matrix(NA,nrow(all.oSPR),2)
    slope.pix <- matrix(NA,nrow(all.oSPR),1)
    for (i in 1:nrow(all.oSPR)){
      print(i)

      nobs <- sum(is.na(all.oSPR[i,])==0 & is.na(all.pSPR[i,])==0)

      if (nobs > 5){
        rmse.pix[i] <- sqrt(mean((all.oSPR[i,]-all.pSPR[i,])^2,na.rm=TRUE))

        corr <- cor(all.oSPR[i,],all.pSPR[i,],use='pairwise.complete.obs')
        corr.test <- cor.test(all.oSPR[i,],all.pSPR[i,])
        cor.pix[i,] <- c(corr,corr.test$p.value)

        slope.pix[i] <- sd(all.oSPR[i,],na.rm=TRUE)/sd(all.pSPR[i,],na.rm=TRUE)
      }
    }

    good.pix[m,1] <- length(which(rmse.pix<7 & cor.pix[,2]<0.05 & slope.pix < 1.2 & slope.pix > 0.8))/nrow(all.oSPR)

    all.oSPR2 <- all.oSPR
    all.pSPR2 <- all.pSPR
    dim(all.oSPR2) <- c(ncell(all.oSPR2),1)
    dim(all.pSPR2) <- c(ncell(all.pSPR2),1)

    #Calculate overlap-wide statistics
    rmse[m,1] <- sqrt(mean((all.oSPR2-all.pSPR2)^2,na.rm=TRUE))
    mbe[m,1] <- mean(all.oSPR2-all.pSPR2,na.rm=TRUE)
    cor[m,1] <- cor(all.oSPR2,all.pSPR2,use='pairwise.complete.obs')

    #Calculate statistics for anomalously early observations
    o_mean <- apply(all.oSPR,1,mean,na.rm=TRUE)
    o_sd <- apply(all.oSPR,1,sd,na.rm=TRUE)

    o_mean_all <- replicate(32, o_mean)
    o_sd_all <- replicate(32, o_sd)
    o_Zscore <- (all.oSPR-o_mean_all)/o_sd_all

    w <- which(o_Zscore < -2)

    rmse_2sd[m,1] <- sqrt(mean((all.oSPR[w]-all.pSPR[w])^2,na.rm=TRUE))
    mbe_2sd[m,1] <- mean(all.oSPR[w]-all.pSPR[w],na.rm=TRUE)
    cor_2sd[m,1] <- cor(all.oSPR[w],all.pSPR[w],use='pairwise.complete.obs')
  }

  save(good.pix,rmse,mbe,cor,rmse_2sd,mbe_2sd,cor_2sd,file = "daymet2landsat_pcam_s50_scene_stats_v4")
})
