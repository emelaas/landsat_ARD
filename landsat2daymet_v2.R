system.time({
  library(ncdf4)
  require(rgdal)
  library(raster)
  library(foreach)
  library(iterators)
  library(doParallel)
  
  args = commandArgs(trailingOnly=T)
  tile_name = args[1]
  
  #tile_name <- 'h20v14'
  
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,sep=''))
  dir.create(paste(getwd(),'/PHENO_1KM/',sep=''))
  
  #Register the parallel backend
  registerDoParallel(16)
  
  deg2rad <- function(deg) {(deg * pi) / (180)}
  
  #Import annual phenology dates
  lc <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',
    tile_name,'/NLCD/nlcd.tif',sep=""))
  w <- c(11,51,52,71,72,73,74,81,82)
  
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
  
  phen <- matrix(NA,ncell(lc),38)
  
  for (j in 1:length(tmp_files)){
    print(j)
    load(tmp_files[j])
    
    phen[(25*5000*(chunk[j]-1)+1):((25*5000)*chunk[j]),] <- all_pheno[,3:40]    
  }
  
  phen[which(getValues(lc) %in% w),] <- NA
  pheno_hdr <- setValues(lc_stack,phen)
  
  
  # Find all Daymet tiles overlapping with ARD tile
  setwd('/projectnb/modislc/data/daymet')
  lon <- raster('daymet_v3_prcp_1981_na.nc4',varname='lon')
  lat <- raster('daymet_v3_prcp_1981_na.nc4',varname='lat')
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
    
    #Convert overlap raster to polygon in order to extract 1km spatial mean phenology dates
    tmp_poly <- rasterToPolygons(cells.crop)
    tmp_poly_proj <- spTransform(tmp_poly,projection(lc))
    
    #Zonal stats for LTM spring/autumn and each inidiv. year
    poly_medians <- foreach(j = 1:length(tmp_poly), .combine = rbind) %dopar% {
      print(j)
      
      tmp_poly_subset <- tmp_poly_proj[j, ]
      poly <- extract(pheno_hdr,tmp_poly_subset)
      poly_mat <- poly[[1]]
      
      if (ncell(poly_mat)>1){
        nobs <- colSums(!is.na(poly_mat))
        poly_mat[,which(nobs<10)] <- NA
        
        poly_median <- round(apply(poly_mat,2,median,na.rm=TRUE))
      } else {
        poly_median <- NA*seq(1,38)
      }
      
    }
    
    obs.SPR <- cbind(cells_vals,poly_medians)
    
    setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
    save(cells.crop,obs.SPR,file = "landsat2daymet")
    
  } else {
    
    #if Landsat scene and Daymet tile extents do not overlap
    obs.SPR <- matrix(NA,ncell(cells.crop),39)
    
    setwd(paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/PHENO_1KM/',sep=""))
    save(cells.crop,obs.SPR,file = "landsat2daymet")
  }
  
})
