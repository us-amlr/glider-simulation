# was 'R_2_rnd_gldrs_Plots&Tables.txt'
# removed '.biom' from plots, fit calculations
# code from 'C:\zot\glider\2019\10oct\oct1\SA_600reps\1_glider'
# calculates fit metrics from 'rndm.means.biom' (glider biomass from
# UFFS.dat) and 'ship.means.biom' (ship). 'biom.gldr' and 'biom.ship' 
# are the mean values from the glider and ship biomasses, respectively.
# copied from jul22, reduced to 2001, 150m SA 
# removed hashes from qntl.vals, etc.
# removed plotting code
# commented-out writes

Fig7 <- function(
  NASC.yrs = c(2001:2009,2011),AMLR.area = 'SA',leg = 1,n.rep = 9,
  gldr.pths = c(1), # one set of glider combinations at a time
  qntl.vals = c(0.97,0.98,0.99,0.999,1,'sum','sd'),
  path.in = paste('tables/',#AMLR.area,'/',
                  sep=''),
  #path.out = paste('tables/',#AMLR.area,'/',
  #                gldr.pths,sep=''), 
  depths = c(150,200,300,400,500,700,1000) # maximum yo depths
  ){
  if(!dir.exists('Fig7'))
    dir.create('Fig7')
  gldr.smpls.m <- gldr.smpls.sd <- array(dim=c(length(depths),length(NASC.yrs)))
  imprecision.table <- bias.table <- array(dim=c(length(depths),length(NASC.yrs)))
  dimnames(imprecision.table) <- dimnames(bias.table) <- list(depths,NASC.yrs)
  for(i.depth in 1:length(depths)){
    yo.count <- read.table(paste('tables/yo_count/',gldr.pths,'_', AMLR.area,'_',
                depths[i.depth],'_yo_count.txt',sep=''))
    ship.sums <- as.matrix(read.table(paste('tables/',
             'shiptable',AMLR.area,'.txt',sep=''),header=TRUE)[,2:8])
    gldr.smpls.yrs <- array(dim=c(n.rep,length(NASC.yrs)))
    dimnames(gldr.smpls.yrs) <- list(1:n.rep,NASC.yrs)
    mae.table <- rmse.table <-array(dim=c(length(NASC.yrs)))
    dimnames(mae.table) <- dimnames(rmse.table) <- list(NASC.yrs)
    for(iyr in 1:length(NASC.yrs)){
      gldr.smpls <- read.table(paste('tables/yo_sums/',gldr.pths,'_', AMLR.area,
                    '_',NASC.yrs[iyr], '_',depths[i.depth],'_yo_sums.txt',sep=''),header=TRUE)
      gldr.smpls.yrs[,iyr] <- apply(as.matrix(gldr.smpls),2,mean,na.rm=TRUE)
      mae.tmp <-  rmse.tmp <- array(dim=c(n.rep))
      for(i in 1:n.rep){
      if(n.rep > 1){
        mae.tmp[i] <- abs(mean(unlist(gldr.smpls[,i]))-ship.sums[iyr,'sum'])/
                                 length(gldr.smpls[,i])
        rmse.tmp[i] <- (mean(unlist(gldr.smpls[,i]))-ship.sums[iyr,'sum'])^2
        }
      else {
        mae.tmp[i] <- as.numeric(abs(mean(unlist(gldr.smpls[,i]))-ship.sums[iyr,'sum']))
        rmse.tmp[i] <- as.numeric(sqrt((mean(unlist(gldr.smpls[i,iyr])-ship.sums[iyr,'sum'])^2)))
        }
      } # end i
    gldr.smpls.m[i.depth,iyr] <- mean(gldr.smpls.yrs[,iyr])
    gldr.smpls.sd[i.depth,iyr] <- sd(gldr.smpls.yrs[,iyr])
    mae.table[iyr] <- sum(mae.tmp) 
    rmse.table[iyr] <- sqrt(sum(rmse.tmp)/
                                 length(gldr.smpls.yrs[,iyr])) 
    } # end iyr

  mae.table <- cbind(mae.table, yo.count = as.numeric(yo.count))
  rmse.table <- cbind(rmse.table, yo.count = as.numeric(yo.count))
  #dimnames(mae.table) <- dimnames(rmse.table) <- list(NASC.yrs,qntl.vals)
  dimnames(yo.count) <- list('yo.count',NASC.yrs)
  imprecision.table[i.depth,]<- round(gldr.smpls.sd[i.depth,]/ship.sums[,'sum'],3)
  bias.table[i.depth,] <- round((gldr.smpls.m[i.depth,]-ship.sums[,'sum'])/ship.sums[,'sum'],3)

  write.table(round(mae.table,1),paste('Fig7/',gldr.pths,'_', AMLR.area,'_',
              depths[i.depth],'_annual_mae.txt',sep=''))
  write.table(round(rmse.table,1),paste('Fig7/',gldr.pths,'_',AMLR.area,'_',depths[i.depth],
              '_annual_rmse.txt',sep=''))
  write.table(imprecision.table,paste('Fig7/',gldr.pths,'_', AMLR.area,'_',
              depths[i.depth],'imprecision.txt',sep=''))
  write.table(bias.table,paste('Fig7/',gldr.pths,'_', AMLR.area,'_',depths[i.depth],'bias.txt',sep=''))
  #######################
  plt.name <- paste('Fig7/',gldr.pths,'_',depths[i.depth],'_Glider & Ship MAE.pdf',sep='')
  pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,oma=c(2,2,2,0))
  y.lim <- c(min(c(mae.table[,'mae.table'],rmse.table[,'rmse.table']),na.rm=TRUE),
           max(c(mae.table[,'mae.table'],rmse.table[,'rmse.table']),na.rm=TRUE))
  plot(as.integer(NASC.yrs),mae.table[,'mae.table'],type='l',col='red',
    ylim=y.lim,main=paste(depths[i.depth],'m',sep=''),
    ylab = 'MAE (red) & RMSE (blue)',xlab = 'Year')
  lines(as.integer(NASC.yrs),rmse.table[,'rmse.table'],col='blue')
  points(as.integer(NASC.yrs),mae.table[,'mae.table'],col='red',pch=19)
  points(as.integer(NASC.yrs),rmse.table[,'rmse.table'],col='blue',pch=19)
  dev.off()

  #######################
  plt.name <- paste('Fig7/',gldr.pths,'_',depths[i.depth],'_Glider & Ship Sv_old.pdf',sep='')
  pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,oma=c(2,2,2,0))
  y.lim <- c(min(c(ship.sums[,'sum'],
                   gldr.smpls.yrs),na.rm=TRUE),
             max(c(ship.sums[,'sum'],
                   gldr.smpls.yrs),na.rm=TRUE))
  plot(NASC.yrs,ship.sums[,'sum'],type='l',col='red',
    ylim=y.lim,lwd=3,cex.lab=0.9,main=depths[i.depth],
    ylab = 'Sv density',xlab = 'Year')
  for(i.rep in 1:n.rep){
    lines(NASC.yrs,gldr.smpls.yrs[i.rep,],col='gray',lwd=3)
    points(NASC.yrs,ship.sums[,'sum'],col='red',pch=19)
    #points(NASC.yrs,gldr.smpls[i.rep,'mean',],col='gray',pch=19)
    }
  lines(NASC.yrs,gldr.smpls.m[i.depth,],col='blue',lwd=3)
  segments(NASC.yrs,gldr.smpls.m[i.depth,] + gldr.smpls.sd[i.depth,],
           NASC.yrs,gldr.smpls.m[i.depth,] - gldr.smpls.sd[i.depth,],
           col='blue',lwd=3)
  points(NASC.yrs,gldr.smpls.m[i.depth,],
         pch=19,col='blue')
  lines(NASC.yrs,ship.sums[,'sum'],col='red',lwd=3)
  dev.off()

  #######################
  plt.name <- paste('Fig7/',gldr.pths,'_',depths[i.depth],'_Glider & Ship Sv_new.pdf',sep='')
  pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,oma=c(2,2,2,0))
  y.lim <- c(min(c(ship.sums[,'sum'],
                   gldr.smpls.yrs),na.rm=TRUE),
             max(c(ship.sums[,'sum'],
                   gldr.smpls.yrs),na.rm=TRUE))
  plot(NASC.yrs,ship.sums[,'sum'],type='l',col='red',
    ylim=y.lim,lwd=3,cex.lab=0.9,main=depths[i.depth],
    ylab = 'Acoustic density',xlab = 'Year')
  for(i.rep in 1:n.rep){
    lines(NASC.yrs,gldr.smpls.yrs[i.rep,],col='gray',lwd=3)
    points(NASC.yrs,ship.sums[,'sum'],col='red',pch=19)
    #points(NASC.yrs,gldr.smpls[i.rep,'sum',],col='gray',pch=19)
    }
  lines(NASC.yrs,gldr.smpls.m[i.depth,],col='blue',lwd=3)
  log.cv <- exp(sqrt(log(1+(gldr.smpls.sd[i.depth,]/gldr.smpls.m[i.depth,])^2)))
  segments(NASC.yrs,gldr.smpls.m[i.depth,] * log.cv,
           NASC.yrs,gldr.smpls.m[i.depth,] / log.cv,
           col='blue',lwd=3)
  points(NASC.yrs,gldr.smpls.m[i.depth,],
         pch=19,col='blue')
  lines(NASC.yrs,ship.sums[,'sum'],col='red',lwd=3)
  dev.off()

  } # end of i.depth loop



} # end Fig7()

