setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/h21v10/IMG')
d <- dir()

home <- '/projectnb/modislc/projects/landsat_sentinel/ARD/h21v10/IMG'

for (i in 1088:length(d)){
  print(i)
  setwd(paste(home,'/',d[i],sep=''))
  d2 <- dir()
  
  n1 <- as.numeric(substr(d[i],16,23))
  n2 <- as.numeric(substr(d[(i+1)],16,23))
  w <- which(as.numeric(substr(d2,16,23))==n2)
  
  for (j in w){
    file.copy(paste(getwd(),'/',d2[j],sep=''),
      paste(home,'/',d[(i+1)],'/',d2[j],sep=''),overwrite=TRUE)
    unlink(paste(getwd(),'/',d2[j],sep=''))
  }
}