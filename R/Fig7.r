# reads Jan12 data
# from C:\zot\glider\2021\10oct\oct25\oct4_vals\oct22
# was 'R_2_rnd_gldrs_Plots&Tables.txt'
# replaced 'gldr.smpls.sd[i.depth,iyr] <- sd(gldr.smpls.yrs[,'sd',iyr])'
# with 'gldr.smpls.sd[i.depth,iyr] <- sd(gldr.smpls.yrs[,'sum',iyr])'

Fig7 <- function(
  NASC.yrs = c(2001:2009,2011),AMLR.area = 'SA',n.rep = 9,
  n.gldr = c(1), # one set of glider combinations at a time
  qntl.vals = c(0.97,0.98,0.99,0.999,1,'sum','sd'),
  depths = c(150,200,300,400,500,700,1000) # maximum yo depths
  ){
  if(!dir.exists('Fig7'))
    dir.create('Fig7')
  gldr.smpls.m <- gldr.smpls.sd <- array(dim=c(length(depths),length(NASC.yrs)))
  imprecision.table <- bias.table <- array(dim=c(length(depths),length(NASC.yrs)))
  dimnames(imprecision.table) <- dimnames(bias.table) <- list(depths,NASC.yrs)
  for(i.depth in 1:length(depths)){
    yo.count <- read.table(paste('tables/yo_count/',n.gldr,'_', AMLR.area,'_',
                depths[i.depth],'_yo_count.txt',sep=''))
    ship.sums <- as.matrix(read.table(paste('tables/',
             'shiptable',AMLR.area,'.txt',sep=''),header=TRUE)[,2:8])
    gldr.smpls.yrs <- array(dim=c(n.rep,length(qntl.vals),length(NASC.yrs)))
    dimnames(gldr.smpls.yrs) <- list(1:n.rep,c(qntl.vals),NASC.yrs)
    mae.table <- rmse.table <- array(dim=c(length(NASC.yrs),length(qntl.vals)))
    dimnames(mae.table) <- dimnames(rmse.table) <- list(NASC.yrs,c(qntl.vals))
    for(iyr in 1:length(NASC.yrs)){
      gldr.smpls <- read.table(paste('tables/gldrtables/',n.gldr,'_', AMLR.area,
                    '_',NASC.yrs[iyr], '_',depths[i.depth],'m_gldrtable.txt',sep=''),header=TRUE)
      gldr.smpls.yrs[,,iyr] <- as.matrix(gldr.smpls)
      mae.tmp <-  rmse.tmp <- array(dim=c(n.rep,length(qntl.vals)))
      for(i in 1:n.rep){
      if(n.rep > 1){
        mae.tmp[i,] <- abs(unlist(gldr.smpls[i,,iyr])-ship.sums[iyr,])/
                                 nrow(gldr.smpls[,,iyr])
        rmse.tmp[i,] <- (unlist(gldr.smpls[i,,iyr])-ship.sums[iyr,])^2
        }
      else {
        mae.tmp[i,] <- as.numeric(abs(unlist(gldr.smpls[i,,iyr])-ship.sums[iyr,]))
        rmse.tmp[i,] <- as.numeric(sqrt((mean(unlist(gldr.smpls[i,,iyr])-ship.sums[iyr,])^2)))
        }
      } # end i
    gldr.smpls.m[i.depth,iyr] <- mean(gldr.smpls.yrs[,'sum',iyr])
    gldr.smpls.sd[i.depth,iyr] <- sd(gldr.smpls.yrs[,'sum',iyr])
    mae.table[iyr,] <- apply(mae.tmp,2,sum) 
    rmse.table[iyr,] <- sqrt(apply(rmse.tmp,2,sum)/
                                 nrow(gldr.smpls.yrs[,,iyr]))  
    } # end iyr

  mae.table <- cbind(mae.table, yo.count = as.numeric(yo.count))
  rmse.table <- cbind(rmse.table, yo.count = as.numeric(yo.count))
  dimnames(yo.count) <- list('yo.count',NASC.yrs)
  imprecision.table[i.depth,]<- round(gldr.smpls.sd[i.depth,]/ship.sums[,'sum'],3)
  bias.table[i.depth,] <- round((gldr.smpls.m[i.depth,]-ship.sums[,'sum'])/ship.sums[,'sum'],3)

  write.table(round(mae.table,1),paste('Fig7/',n.gldr,'_', AMLR.area,'_',
              depths[i.depth],'_annual_mae.txt',sep=''))
  write.table(round(rmse.table,1),paste('Fig7/',n.gldr,'_',AMLR.area,'_',depths[i.depth],
              '_annual_rmse.txt',sep=''))
  write.table(imprecision.table,paste('Fig7/',n.gldr,'_', AMLR.area,'_',
              depths[i.depth],'imprecision.txt',sep=''))
  write.table(bias.table,paste('Fig7/',n.gldr,'_', AMLR.area,'_',depths[i.depth],'bias.txt',sep=''))
  #######################
  plt.name <- paste('Fig7/',n.gldr,'_',AMLR.area,'_',depths[i.depth],
                    '_Glider & Ship MAE.pdf',sep='')
  pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,oma=c(2,2,2,0))
  y.lim <- c(min(c(mae.table[,'sum'],rmse.table[,'sum']),na.rm=TRUE),
           max(c(mae.table[,'sum'],rmse.table[,'sum']),na.rm=TRUE))
  plot(as.integer(NASC.yrs),mae.table[,'sum'],type='l',col='red',
    ylim=y.lim,main=paste(depths[i.depth],'m',sep=''),
    ylab = 'MAE (red) & RMSE (blue)',xlab = 'Year')
  lines(as.integer(NASC.yrs),rmse.table[,'sum'],col='blue')
  points(as.integer(NASC.yrs),mae.table[,'sum'],col='red',pch=19)
  points(as.integer(NASC.yrs),rmse.table[,'sum'],col='blue',pch=19)
  dev.off()


  #######################
  plt.name <- paste('Fig7/',n.gldr,'_',AMLR.area,'_',depths[i.depth],
                    '_Glider & Ship Sv.pdf',sep='')
  pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,oma=c(2,2,2,0))
  y.lim <- c(min(c(ship.sums[,'sum'],
                   gldr.smpls.yrs[,6,]),na.rm=TRUE),
             max(c(ship.sums[,'sum'],
                   gldr.smpls.yrs[,6,]),na.rm=TRUE))
  plot(NASC.yrs,ship.sums[,'sum'],type='l',col='red',
    ylim=y.lim,lwd=3,cex.lab=0.9,main=depths[i.depth],
    ylab = 'Acoustic density',xlab = 'Year')
  for(i.rep in 1:n.rep){
    lines(NASC.yrs,gldr.smpls.yrs[i.rep,'sum',],col='gray',lwd=3)
    points(NASC.yrs,ship.sums[,'sum'],col='red',pch=19)
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

