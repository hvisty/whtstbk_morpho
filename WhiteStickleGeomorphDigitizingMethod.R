library(geomorph)
setwd("/Users/Hannah/Documents/Schluter Lab Work/FishPictures")

filelist<-list.files("/Users/Hannah/Documents/Schluter Lab Work/FishPictures")
filelist

digitize2d(filelist, nlandmarks=8, scale=1, tpsfile="white_stickleback_morphometry.tps", verbose=F)

