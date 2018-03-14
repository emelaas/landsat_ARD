library(stringr)
library(RCurl)

#args = commandArgs(trailingOnly=T)
#tile_name <- args[1]

tile_name <- 'h18v13'

setwd('/projectnb/modislc/projects/landsat_sentinel/ARD/')
dir.create(paste(getwd(),'/',tile_name,sep=''))

url <- paste("https://edclpdsftp.cr.usgs.gov/downloads/collections/l2-ard-tiles/",
  tile_name,"/",sep='')
html <- paste(readLines(url), collapse="\n")
matched <- str_match_all(html, "<a href=\"(.*?)\"")

tmp <- matched[[1]]

filenames <- tmp[,2]
ta <- grep("TA.tar",filenames)
sr <- grep("SR.tar",filenames)

for (i in 1:length(sr)){
  print(i)
  download.file(paste(url,filenames[sr[i]],sep=''),
    destfile=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',
    filenames[sr[i]],sep=''),method="wget",quiet=TRUE)
  download.file(paste(url,filenames[ta[i]],sep=''),
    destfile=paste('/projectnb/modislc/projects/landsat_sentinel/ARD/',tile_name,'/',
    filenames[ta[i]],sep=''),method="wget",quiet=TRUE)
}
