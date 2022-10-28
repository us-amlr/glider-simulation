# this function combines annual tables produced by 'Fig7.r' into
# rmse.means, mae.means, bias.table and imprecision.table
# and produces Figs 7 and 8 contour plots of these values

Fig7.8contours <- function(
  NASC.yrs =c(2001:2009,2011),
  AMLR.area = 'SA',
  n.reps = 9,
  n.gldr = c(1,2,3,4,5),
  depths = c(150,200,300,400,500,700,1000)){

  library(viridis)
  gldr.smpls.m <- gldr.smpls.sd <- array(dim=c(length(depths),length(NASC.yrs)))
  mae.means <- rmse.means <- sd.table <- bias.table <- imprecision.table <-
    array(dim=c(length(n.gldr),length(depths)))
  dimnames(mae.means) <- dimnames(rmse.means) <- dimnames(sd.table) <-
    dimnames(bias.table) <- dimnames(imprecision.table) <- list(n.gldr,depths)

  sd.table <- array(dim=c(length(depths),length(n.gldr)))
  dimnames(sd.table) <- list(depths,n.gldr)
  for(i.gldr in 1: length(n.gldr))
    for(i.depth in 1:length(depths)){
      path.bias <- paste('Fig7/',n.gldr[i.gldr],'_',AMLR.area,'_',depths[i.depth],sep='')

      mae.table <- read.table(paste('Fig7/',n.gldr[i.gldr],'_', AMLR.area,
            '_',depths[i.depth],'_annual_mae.txt',sep=''))
      rmse.table <- read.table(paste('Fig7/',n.gldr[i.gldr],'_', AMLR.area,
            '_',depths[i.depth],'_annual_rmse.txt',sep='')) 
      mae.means[i.gldr,i.depth] <- mean(mae.table[,'sum'])
      rmse.means[i.gldr,i.depth] <- mean(rmse.table[,'sum'])

      bias.table[i.gldr,] <-
            apply(read.table(paste(path.bias,'bias.txt',sep=''),header=TRUE),1,mean)
      imprecision.table[i.gldr,] <- apply(read.table(paste(path.bias,'imprecision.txt',sep=''),
	                        header=TRUE),1,mean)
    } # end i.gldr, i.depth

  if(length(n.gldr)>=3){ # need 3 or more combined glider samples/replicate for contouring
    par(cex=2,cex.lab=1.8)   #,mar=c(5,6,4,1)+.05)
    filled.contour(y=depths,x=1:length(n.gldr),z=imprecision.table,
        plot.axes = {
          axis(1)
          axis(2)
          contour(y=depths,x=1:length(n.gldr),z=imprecision.table, add=TRUE,labcex=2)
          },
        color.palette = viridis,
        plot.title = title(main = 'Imprecision\n',cex.lab=1.8,ylab = 'Max Depth (m)',xlab='N gliders'),
        key.title=title(main='(m)'),
        nlevels=12)
	points(x=rep(1:5,length(depths)),y=rep(depths,5))
    #dev.new()
    #jpeg(paste('Fig8bias_',AMLR.area,'.jpeg',sep=''))
    par(cex=1.5,cex.lab=1.8)
    filled.contour(y=depths,x=1:length(n.gldr),z=bias.table,
        plot.axes = {
          axis(1)
          axis(2)
          contour(y=depths,x=1:length(n.gldr),z=bias.table, add=TRUE,labcex=2)
          },
        color.palette = viridis,
        plot.title = title(main = 'Bias\n',cex.lab=1.8,ylab = 'Max Depth (m)',xlab='N gliders'),
        key.title=title(main='(m)'),
        nlevels=12)
	points(x=rep(1:5,length(depths)),y=rep(depths,5))
    #dev.new()
    #jpeg(paste('Fig7RMSE_',AMLR.area,'.jpeg',sep=''))
    par(cex=1.5,cex.lab=1.8)
    filled.contour(y=depths,x=1:length(n.gldr),z=rmse.means,
        plot.axes = {
          axis(1)
          axis(2)
          contour(y=depths,x=1:length(n.gldr),z=rmse.means, add=TRUE,labcex=2)
          },
        color.palette = viridis,
        plot.title = title(main = 'RMSE\n',cex.lab=1.8,ylab = 'Max Depth (m)',xlab='N gliders'),
        key.title=title(main='(m)'),
        xlab='N gliders',ylab = 'Max Depth (m)',nlevels=12)
	points(x=rep(1:5,length(depths)),y=rep(depths,5))
    #dev.new()
    #jpeg(paste('Fig7MAE_',AMLR.area,'.jpeg',sep=''))
    par(cex=1.5,cex.lab=1.8)
    filled.contour(y=depths,x=1:length(n.gldr),z=mae.means,
        plot.axes = {
          axis(1)
          axis(2)
          contour(y=depths,x=1:length(n.gldr),z=mae.means, add=TRUE,labcex=2)
          },
        color.palette = viridis,
        plot.title = title(main = 'MAE\n',cex.lab=1.8,ylab = 'Max Depth (m)',xlab='N gliders'),
        key.title=title(main='(m)'),
        xlab='N gliders',ylab = 'Max Depth (m)',nlevels=12)
	points(x=rep(1:5,length(depths)),y=rep(depths,5))
        } 
  } # end Fig7_contours()
