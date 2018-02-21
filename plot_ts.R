i <- which.max(all_pheno[,4])
evi_ts <- as.matrix(dt.evi[,..i]/10000)

x11()
par(mfrow=c(6,6),mar=c(0,0,0,0))
for (y in 1982:2017){
  plot(doy,evi_ts,pch=16,col='gray')
  points(doy[yr==y],evi_ts[yr==y],pch=21,col='black',bg='orange',cex=1.3)
  abline(v=all_pheno[i,3],lty=2)
}
