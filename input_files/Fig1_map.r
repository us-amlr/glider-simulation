library(maps)
library(mapdata)
library(splancs)
library(rgdal)

Fig1 <- function(){
  SA.lines <- read.csv('input_files/AMLR_SA_transects.csv')
  WA.lines <- read.csv('input_files/AMLR_WA_transects.csv')
  SA.x <- SA.lines$Long
  SA.y <- SA.lines$Lat
  SA.labs <- SA.lines$Station
  WA.x <- WA.lines$Long
  WA.y <- WA.lines$Lat
  WA.labs <- WA.lines$Station

  #jpeg('Fig1_bw.jpeg')
  map('worldHires',xlim=c(-73,-51), ylim=c(-65,-54.5), fill=T,col='gray',ylab='Latitude')
  map.axes()
  text(-71,-55.5,substitute(paste(bold('Chile'))),srt=-27,cex=1.3)
  text(-64.5,-55.7,substitute(paste(bold('Argentina'))),cex=1.3)  
  text(-62,-61.35,substitute(paste(bold('Cape Shirreff'))),srt=45,lwd=3) 
  text(-54,-62.1,substitute(paste(bold('Bransfield Strait'))),lwd=3) 
  lines(cbind(x=WA.x,y=WA.y),lwd=2)
  lines(cbind(x=SA.x,y=SA.y),lwd=2)
  #dev.off()
} # end Fig1