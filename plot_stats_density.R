below5 <- raster('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs/below5_sum.tif')
below5_vals <- getValues(below5)
w1 <- which(below5_vals<=2)
w2 <- which(below5_vals>2 & below5_vals<6)
w3 <- which(below5_vals>=6)

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("cor*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals1 <- s_vals[w1,]
s_vals2 <- s_vals[w2,]
s_vals3 <- s_vals[w3,]
colnames(s_vals) <- c('Chill','Photo','SW Mar 17','SW Jan1 (Mean AGDD)')
dens1 <- apply(s_vals1, 2, density, na.rm=TRUE)
dens2 <- apply(s_vals2, 2, density, na.rm=TRUE)
dens3 <- apply(s_vals3, 2, density, na.rm=TRUE)

#x11(h=4,w=11)
pdf(h=4,w=11,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/cor_density_all.pdf')
par(mfrow=c(1,3))
plot(NA, xlim=c(0,1), ylim=range(sapply(dens1, "[", "y")),xlab='Correlation Coefficient',ylab='',
  main=expression('0-2 months with T'[avg]*' < 5'~degree*'C'))
mapply(lines, dens1, col=1:length(dens1),lwd=4)
plot(NA, xlim=c(0,1), ylim=range(sapply(dens2, "[", "y")),xlab='Correlation Coefficient',ylab='',
  main=expression('3-5 months with T'[avg]*' < 5'~degree*'C'))
mapply(lines, dens2, col=1:length(dens2),lwd=4)
plot(NA, xlim=c(0,1), ylim=range(sapply(dens3, "[", "y")),xlab='Correlation Coefficient',ylab='',
  main=expression('6-8 months with T'[avg]*' < 5'~degree*'C'))
mapply(lines, dens3, col=1:length(dens3),lwd=4)

legend("topleft", legend=colnames(s_vals), fill=1:length(dens3),bty='n')
dev.off()


setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("RMSE*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals1 <- s_vals[w1,]
s_vals2 <- s_vals[w2,]
s_vals3 <- s_vals[w3,]
colnames(s_vals) <- c('Chill','Photo','SW Mar 17','SW Jan1 (Mean AGDD)')
dens1 <- apply(s_vals1, 2, density, na.rm=TRUE)
dens2 <- apply(s_vals2, 2, density, na.rm=TRUE)
dens3 <- apply(s_vals3, 2, density, na.rm=TRUE)

#x11(h=4,w=11)
pdf(h=4,w=11,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/rmse_density_all.pdf')
par(mfrow=c(1,3))
plot(NA, xlim=c(2,14), ylim=range(sapply(dens1, "[", "y")),xlab='RMSE',ylab='',
  main=expression('0-2 months with T'[avg]*' < 5'~degree*'C'))
mapply(lines, dens1, col=1:length(dens1),lwd=4)
plot(NA, xlim=c(2,14), ylim=range(sapply(dens2, "[", "y")),xlab='RMSE',ylab='',
  main=expression('3-5 months with T'[avg]*' < 5'~degree*'C'))
mapply(lines, dens2, col=1:length(dens2),lwd=4)
plot(NA, xlim=c(2,14), ylim=range(sapply(dens3, "[", "y")),xlab='RMSE',ylab='',
  main=expression('6-8 months with T'[avg]*' < 5'~degree*'C'))
mapply(lines, dens3, col=1:length(dens3),lwd=4)

legend("topright", legend=colnames(s_vals), fill=1:length(dens3),bty='n')
dev.off()


setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/tifs')
in_dirs_tile <- list.files(path=getwd(),
  pattern=glob2rx("slope*tif"),full.names=T,include.dirs=T,recursive=TRUE)
s <- stack(in_dirs_tile)
s_vals <- getValues(s)
s_vals1 <- s_vals[w1,]
s_vals2 <- s_vals[w2,]
s_vals3 <- s_vals[w3,]
colnames(s_vals) <- c('Chill','Photo','SW Mar 17','SW Jan1 (Mean AGDD)')
dens1 <- apply(s_vals1, 2, density, na.rm=TRUE)
dens2 <- apply(s_vals2, 2, density, na.rm=TRUE)
dens3 <- apply(s_vals3, 2, density, na.rm=TRUE)

#x11(h=4,w=11)
pdf(h=4,w=11,'/projectnb/modislc/projects/landsat_sentinel/ARD/figures/NCC/slope_density_all.pdf')
par(mfrow=c(1,3))
plot(NA, xlim=c(0,2), ylim=range(sapply(dens1, "[", "y")),xlab='RMA Slope',ylab='',
  main=expression('0-2 months with T'[avg]*' < 5'~degree*'C'))
mapply(lines, dens1, col=1:length(dens1),lwd=4)
abline(v=1,lty=2)
plot(NA, xlim=c(0,2), ylim=range(sapply(dens2, "[", "y")),xlab='RMA Slope',ylab='',
  main=expression('3-5 months with T'[avg]*' < 5'~degree*'C'))
mapply(lines, dens2, col=1:length(dens2),lwd=4)
abline(v=1,lty=2)
plot(NA, xlim=c(0,2), ylim=range(sapply(dens3, "[", "y")),xlab='RMA Slope',ylab='',
  main=expression('6-8 months with T'[avg]*' < 5'~degree*'C'))
mapply(lines, dens3, col=1:length(dens3),lwd=4)
abline(v=1,lty=2)

legend("topleft", legend=colnames(s_vals), fill=1:length(dens3),bty='n')
dev.off()


