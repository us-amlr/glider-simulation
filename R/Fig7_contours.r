# from 'C:\zot\glider\2021\1jan\jan12\R_5\SA\R_gldr_contour_SA.txt'
# from 'Users/noaa/glider/2020/4apr/apr15_fit_contours/fit_contours_SA/'
# from '/volumes/DKinzey_Backup/glider/4apr/'
# this file combines annual tables into
# rmse.means, mae.means, bias.table and imprecision.table
# and produces Fig 5 contour plots of these values

Fig7_contours <- function(NASC.yrs =c(2001:2009,2011),
AMLR.area = 'SA',
n.reps = 9,
quant.mos = c(0.97,0.98,0.99,0.999,1,'mean','sd'),
n.gldr.pth = c(1,2,3,4,5),
yo.depths = c(150,200,300,400,500,700,1000)){

gldr.smpls.m <- gldr.smpls.sd <- array(dim=c(length(yo.depths),length(NASC.yrs)))
mae.means <- rmse.means <- sd.table <- bias.table <- imprecision.table <-
  array(dim=c(length(n.gldr.pth),length(yo.depths)))
dimnames(mae.means) <- dimnames(rmse.means) <- dimnames(sd.table) <-
  dimnames(bias.table) <- dimnames(imprecision.table) <- list(n.gldr.pth,yo.depths)

sd.table <- array(dim=c(length(yo.depths),length(n.gldr.pth)))
dimnames(sd.table) <- list(yo.depths,n.gldr.pth)
for(i.gldr in 1: length(n.gldr.pth))
  for(i.depth in 1:length(yo.depths)){
    path.bias <- paste('Fig7/',n.gldr.pth[i.gldr],'_',AMLR.area,'_',yo.depths[i.depth],sep='')

    mae.table <- read.table(paste('Fig7/',n.gldr.pth[i.gldr],'_', AMLR.area,
            '_',yo.depths[i.depth],'_annual_mae.txt',sep=''))
    rmse.table <- read.table(paste('Fig7/',n.gldr.pth[i.gldr],'_', AMLR.area,
            '_',yo.depths[i.depth],'_annual_rmse.txt',sep='')) 
    mae.means[i.gldr,i.depth] <- mean(mae.table[,'mae.table'])
    rmse.means[i.gldr,i.depth] <- mean(rmse.table[,'rmse.table'])

  bias.table[i.gldr,] <-
            apply(read.table(paste(path.bias,'bias.txt',sep=''),header=TRUE),1,mean)
  imprecision.table[i.gldr,] <- apply(read.table(paste(path.bias,'imprecision.txt',sep=''),
	                        header=TRUE),1,mean)
  } # end i.gldr, i.depth

if(length(n.gldr.pth)>3){ # need 3 or more glider combinations/replicate for contouring
  par(cex=1.5,cex.lab=1.5)
  contour(y=yo.depths,x=1:length(n.gldr.pth),z=imprecision.table,main = 'Imprecision',
        xlab='N gliders',ylab = 'Max Yo Depth (m)',nlevels=12,labcex=1.5)
	points(x=rep(1:5,length(yo.depths)),y=rep(yo.depths,5))
  #dev.new()
  par(cex=1.5,cex.lab=1.5)
  contour(y=yo.depths,x=1:length(n.gldr.pth),z=bias.table,main = 'Bias',
        xlab='N gliders',ylab = 'Max Yo Depth (m)',nlevels=12,labcex=1.5)
	points(x=rep(1:5,length(yo.depths)),y=rep(yo.depths,5))
  #dev.new()
  par(cex=1.5,cex.lab=1.5)
  contour(y=yo.depths,x=1:length(n.gldr.pth),z=rmse.means,main = 'RMSE',
        xlab='N gliders',ylab = 'Max Yo Depth (m)',nlevels=12,labcex=1.5)
	points(x=rep(1:5,length(yo.depths)),y=rep(yo.depths,5))
  #dev.new()
  par(cex=1.5,cex.lab=1.5)
  contour(y=yo.depths,x=1:length(n.gldr.pth),z=mae.means,main = 'MAE',
        xlab='N gliders',ylab = 'Max Yo Depth (m)',nlevels=12,labcex=1.5)
	points(x=rep(1:5,length(yo.depths)),y=rep(yo.depths,5))
  } # end length(n.gldr.pth)>3
} # end Fig7_contours()
