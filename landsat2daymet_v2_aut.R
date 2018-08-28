system.time({
  library(ncdf4)
  require(rgdal)
  library(raster)
  library(foreach)
  library(iterators)
  library(doParallel)

  args = commandArgs(trailingOnly=T)
  tile_name = args[1]
  #tile_name <- 'h19v05'

  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,sep=''))
  dir.create(paste(getwd(),'/PHENO_1KM/',sep=''))

  #Register the parallel backend
  registerDoParallel(16)

  deg2rad <- function(deg) {(deg * pi) / (180)}

  #Import annual phenology dates
  lc <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/NLCD/nlcd.tif',sep=""))
  w <- c(0,11,21,22,23,24,31,51,52,71,72,73,74,81,82,90,95)

  lc_vals <- getValues(lc)
  lc_vals[lc_vals>=-9999] <- NA
  lc_NA <- setValues(lc,lc_vals)
  lc_stack <- stack(mget(rep("lc_NA",38)))

  ##### Calculate zonal mean of phenology across Daymet polygon map #####

  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/IMG',sep=''))
  in_dirs <- list.files(path=getwd(),pattern=glob2rx("L*"),
    full.names=T,include.dirs=T)

  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO',sep=''))
  tmp_files <- list.files(pattern = "evi2_phenology_topocorr*", recursive = TRUE,full.names = TRUE)
  chunk <- unlist(lapply(tmp_files,
    function(x) na.omit(as.numeric(unlist(strsplit(unlist(x), "[^0-9]+"))))[2]))

  phen_aut <- matrix(NA,ncell(lc),38)
  
  rsmooth <- matrix(NA,ncell(lc),1)

  for (j in 1:length(tmp_files)){
    print(j)
    load(tmp_files[j])

    phen_aut[(25*5000*(chunk[j]-1)+1):((25*5000)*chunk[j]),] <- all_pheno[,c(3:4,41:76)]
    rsmooth[(25*5000*(chunk[j]-1)+1):((25*5000)*chunk[j]),] <- all_pheno[,2]
  }

  phen_aut[which(getValues(lc) %in% w),] <- NA
  phen_aut[which(rsmooth<0.9),] <- NA
  phen_aut_hdr <- setValues(lc_stack,phen_aut)

  writeRaster(phen_aut_hdr[[1]],filename=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/MAPS/aut_mask',sep=""),format='GTiff',overwrite=TRUE)

  # Find all Daymet tiles overlapping with ARD tile
  setwd('/projectnb/modislc/data/climate/daymet')
  tmax <- raster('daymet_v3_tmax_1981_na.nc4',varname='tmax')
  cells <- setValues(tmax,seq(1,ncell(tmax)))

  #Native Daymet projection
  daymet_crs <- CRS("+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +a=6378137 +rf=298.257223563 +lat_1=25 +lat_2=60")

  #Import overlap shapefile
  o <- readOGR(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/SHP/',sep=""),tile_name)

  #Reproject overlap shapefile
  o_reproj <- spTransform(o,daymet_crs)

  cells.crop <- crop(cells,extent(o_reproj))
  cells_vals = getValues(cells.crop)

  #Convert overlap raster to polygon in order to extract 1km spatial mean phenology dates
  tmp_poly <- rasterToPolygons(cells.crop)
  tmp_poly_proj <- spTransform(tmp_poly,projection(lc))

  #Zonal stats for LTM spring/autumn and each inidiv. year
  poly_medians <- foreach(j = 1:length(tmp_poly), .combine = rbind) %dopar% {
    print(j)

    tmp_poly_subset <- tmp_poly_proj[j, ]
    poly <- extract(phen_aut_hdr,tmp_poly_subset)
    poly_mat <- poly[[1]]

    if (ncell(poly_mat)>1){
      nobs <- colSums(!is.na(poly_mat))
      poly_mat[,which(nobs<50)] <- NA

      poly_median <- round(apply(poly_mat,2,median,na.rm=TRUE))
    } else {
      poly_median <- NA*seq(1,38)
    }

  }

  obs.AUT <- cbind(cells_vals,poly_medians)

  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
  save(cells.crop,obs.AUT,file = "landsat2daymet_aut")

})
