library(raster)
library(rgdal)
library(gdalUtils)
library(foreach)
library(iterators)
library(doParallel)
library(ncdf4)
require(RColorBrewer)
require(classInt)

#Register the parallel backend
registerDoParallel(8)

r_poly <- readOGR('/projectnb/modislc/users/emelaas/scratch32/FIA_species_map/RDS-2013-0013_Data/Data/SHPfiles/','easternUS25k')
r_reproj <- spTransform(r_poly,
  CRS('+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'))

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/')
in_dirs <- list.files(path=getwd(),
  pattern=glob2rx('cor*tif'),full.names=T,include.dirs=T,recursive=T)

#Get list of individual bands
cor <- stack(in_dirs)

cor_hex_all <- foreach(i = 1:length(r_reproj), .combine = rbind) %dopar% {
  print(i)    
  
  cor_hex_mean <- c(r_reproj$proj_ID25k[i],matrix(NA,1,length(in_dirs)))
  
  r_extr <- extract(cor,r_reproj[i,])
  if (sapply(r_extr,is.null)==1){
    
  } else if (length(r_extr[[1]][,1])==1){
    
  } else {
    
    cor_hex_mean <- c(r_reproj$proj_ID25k[i],t(as.matrix(colMeans(r_extr[[1]],na.rm=TRUE))))
  }
  
  cor_hex <- cor_hex_mean
}

# Save model calibration results
save(cor_hex_all,file='/projectnb/modislc/projects/landsat_sentinel/ARD/R_files/FIA_hex_cor_mean')

load(file='/projectnb/modislc/projects/landsat_sentinel/ARD/R_files/FIA_hex_cor_mean')

spp_grps <- read.table(file='/projectnb/modislc/users/emelaas/scratch32/FIA_species_map/east_hex_groups.csv',
  sep=',')

all_spp_mat <- matrix(NA,12548,323)
setwd('/projectnb/modislc/users/emelaas/scratch32/FIA_species_map/RDS-2013-0013_Data/Data/RasterMapCIs')
in_dirs <- list.files(path=getwd(),
  pattern=glob2rx('*25k*'),full.names=T,include.dirs=T,recursive=T)
for (i in 1:length(in_dirs)){
  print(i)
  tmp <- read.table(file=in_dirs[i],sep=',',header=TRUE)
  
  matches <- regmatches(in_dirs[i], gregexpr("[[:digit:]]+", in_dirs[i]))
  spp_num <- as.numeric(unlist(matches))[4]
  
  all_spp_mat[,which(spp_grps[,1]==spp_num)] <- tmp[,2]
}

east_spp_mat <- all_spp_mat[r_reproj$proj_ID25k,]
east_tot <- replicate(323,rowSums(east_spp_mat))
east_spp_perc <- east_spp_mat/east_tot
east_spp_perc[is.na(east_spp_perc)==1] <- 0

east_grp_perc <- matrix(NA,6215,5)
for (i in 1:5){
  east_grp_perc[,i] <- rowSums(east_spp_perc[,spp_grps[,3]==i])
}

grp_max <- apply(east_grp_perc,1,which.max)
cor.group <- matrix(NA,5,4)
cor.error.group <- matrix(NA,5,4)
for (i in 1:4){
  cor.group[,i] <- aggregate(cor_hex_all[,(i+1)],by=list(grp_max),FUN=mean,na.rm=TRUE)[,2]
  cor.error.group[,i] <- aggregate(cor_hex_all[,(i+1)],by=list(grp_max),FUN=sd,na.rm=TRUE)[,2]
}

models <- c('Chill','Photo','SW Mar 17','SW Jan 1 (Opt)')
colors <- c('royalblue','springgreen','orange','yellow')

#x11(h=5,w=8)
pdf(h=5,w=8,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/cor_by_ftype.pdf')
barCenters <- barplot(t(cor.group),beside=TRUE,
  names=c('OH','OG','EAC','MBB','AB'),ylim=c(0,1),
  col=colors, ylab='Correlation Coefficient',xlab='Major Forest Type')
segments(barCenters, t(cor.group)-t(cor.error.group), barCenters,
  t(cor.group)+t(cor.error.group), lwd=1.5)
legend('topright',legend=models,
  cex=0.8,fill=colors,
  ncol=2)
dev.off()


