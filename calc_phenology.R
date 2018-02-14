library(rgdal)
library(raster)
library(foreach)
library(iterators)
library(doParallel)
library(zoo)
library(RColorBrewer)
require(data.table)
library(lubridate)

source(file='/usr3/graduate/emelaas/Code/R/landsat_sentinel/v1_3/MCD12Q2C6_functions_1.R')

args = commandArgs(trailingOnly=T)
chunk = as.numeric(args[1])

# chunk <- 1

#Register the parallel backend
registerDoParallel(16)

comb <- function(x, ...) {
  lapply(seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

pheno_pars <- list(
  LandsatFillQuant=0.05,
  LandsatXmin=0,
  LandsatSpikeThresh=2,
  LandsatMinResid=0.1,
  LandsatFillDOY=NULL,
  LandsatDoAnnual=T,
  LandsatPadHeadTail=T,
  min_peak_to_peak_distance=50,
  min_peak_quantile=0.2,
  max_seg_length=200,
  min_seg_amplitude=0.00,
  agg_amp_frac=0.15,
  gup_threshes=c(0.1,0.5,0.9),
  gdown_threshes=c(0.9,0.5,0.1),
  spline_spar=0.5
)
#---------------------------------------------------------------------

for (k in 2:2){
  #Find all images for each site
  setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/h22v15/')
  in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx("evi2.tif"),full.names=T,include.dirs=T,recursive=TRUE)

  r <- stack(in_dirs_tile)
  c <- crop(r,extent(r, 25*(chunk-1)+1, 25*(chunk), 1, 5000))
  vi_vals <- getValues(c)
  vi_vals[vi_vals<=0] <- NA

  #Generate index of ancillary info for all images
  doy <- yday(ymd(substring(in_dirs_tile,72,79)))
  yr <- as.numeric(substring(in_dirs_tile,72,75))
  sat <- substring(in_dirs_tile,60,60)

  index <- data.frame(sat,yr,doy,t(vi_vals))
  colnames(index) <- c('sat','yr','doy','vi')
  index <- index[order(yr,doy),]
  index$time <- strptime(paste(index$yr, index$doy), format="%Y %j")

  #index2 <- index2[-which(index2$sat=='S'),]

  vi_vals <- index[,-c(1,2,3,ncol(index))]
  if (ncol(vi_vals)<10000){
    block_width <- ncol(vi_vals)-1
  } else {
    block_width <- 10000
  }
  nblocks <- ncol(vi_vals)%/%block_width
  bs_start <- seq(1,nblocks*block_width+1,block_width)
  bs_end <- c(seq(block_width,nblocks*block_width,block_width),ncol(vi_vals))

  dt.evi <- data.table(vi_vals,keep.rownames=FALSE)
  info <- as.matrix(index[,1:3])

  system.time({
      pheno_mat_all <- foreach(i = 1:ncol(vi_vals), .combine='comb', .multicombine=TRUE,
        .init=list(list(), list(), list(), list(), list(), list())) %dopar% {
        if (i%%10000==0) print(i)

        all <- list(as.numeric(matrix(NA,14,1)),as.numeric(matrix(NA,14,1)),
          as.numeric(matrix(NA,14,1)),as.numeric(matrix(NA,14,1)),as.numeric(matrix(NA,14,1)))

        #Generate index table
        index2 <- cbind(info,dt.evi[,..i]/10000)
        time <- as.Date(substr(strptime(paste(info[,1], info[,2]), format="%Y %j"),1,10))
        daynum <- as.numeric(time-time[1])
        index3 <- index2[-which(is.na(index2[,4])==1),]
        time <- time[-which(is.na(index2[,4]))]

        pheno_metrics <- try(DoPhenologyLandsat(as.matrix(index3[,4]),time,pheno_pars))

        if (nrow(index3)>0 & is.na(pheno_metrics)==0){
          num_obs <- nrow(index3) #Calculate number of cloud-free observations
          w_gup <- which(pheno_metrics[[8]]-pheno_metrics[[7]]>0.1) #which green up shoulders have amplitude > 0.10?
          w_gdown <- which(pheno_metrics[[10]]-pheno_metrics[[9]]>0.1) #which green down shoulders?

          num_cyc_gup <- length(w_gup) #total number of green ups
          num_cyc_gdown <- length(w_gdown) #total number of green downs

          if (min(c(num_cyc_gup,num_cyc_gdown))>0){
            all_yrs <- as.numeric(substr(pheno_metrics[[1]],1,4))
            w_gup_yrs <- all_yrs[w_gup]-2012
            w_gdown_yrs <- all_yrs[w_gdown]-2012


            # spr1, spr2, spr3, seg_min_gup, seg_max_gup, gup_rsq, gup_missing, gup_longest
            tmp <- unlist(pheno_metrics)
            for (j in 1:length(w_gup)){
              a <- tmp[seq(w_gup[j],length(tmp),length(all_yrs))] #green up metrics
              all[[w_gup_yrs[j]]][c(1:8)] <- a[c(1,2,3,7,8,12,13,14)]
            }

            # aut1, aut2, aut3, gdown_rsq, gdown_missing, gdown_longest
            for (j in 1:length(w_gdown)){
              b <- tmp[seq(w_gdown[j],length(tmp),length(all_yrs))] #green down metrics
              all[[w_gdown_yrs[j]]][c(9:14)] <- b[c(4,5,6,15,16,17)]
            }
          } else {
            num_obs <- nrow(index3)
          }

        } else {
          num_obs <- nrow(index3)
          num_cyc_gup <- 0
          num_cyc_gdown <- 0
        }

        list(c(num_obs,num_cyc_gup,num_cyc_gdown),all[[1]],all[[2]],all[[3]],all[[4]],all[[5]])
      }

  })

  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',site$code,sep=""))
  save(pheno_mat_all,file=paste(site$code,'_','evi2_corr_phenology_',chunk,sep=""))
}
