# this function combines annual tables into
# rmse.means, mae.means, bias.table and imprecision.table
# and produces Fig 7 contour plots of these values

Fig7_contours <- function(
  NASC.yrs =c(2001:2009,2011),
  AMLR.area = 'SA',
  n.reps = 9,
  n.gldr = c(1,2,3,4,5),
  depths = c(150,200,300,400,500,700,1000)){

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
    par(cex=1.5,cex.lab=1.5)
    contour(y=depths,x=1:length(n.gldr),z=imprecision.table,main = 'Imprecision',
        xlab='N gliders',ylab = 'Max Yo Depth (m)',nlevels=12,labcex=1.5)
	points(x=rep(1:5,length(depths)),y=rep(depths,5))
    #dev.new()
    par(cex=1.5,cex.lab=1.5)
    contour(y=depths,x=1:length(n.gldr),z=bias.table,main = 'Bias',
        xlab='N gliders',ylab = 'Max Yo Depth (m)',nlevels=12,labcex=1.5)
	points(x=rep(1:5,length(depths)),y=rep(depths,5))
    #dev.new()
    par(cex=1.5,cex.lab=1.5)
    contour(y=depths,x=1:length(n.gldr),z=rmse.means,main = 'RMSE',
        xlab='N gliders',ylab = 'Max Yo Depth (m)',nlevels=12,labcex=1.5)
	points(x=rep(1:5,length(depths)),y=rep(depths,5))
    #dev.new()
    par(cex=1.5,cex.lab=1.5)
    contour(y=depths,x=1:length(n.gldr),z=mae.means,main = 'MAE',
        xlab='N gliders',ylab = 'Max Yo Depth (m)',nlevels=12,labcex=1.5)
	points(x=rep(1:5,length(depths)),y=rep(depths,5))
    } # end length(n.gldr)>3
} # end Fig7_contours()
