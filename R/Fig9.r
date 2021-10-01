# from 'C:\zot\glider\2021\2feb\feb2_depth_contours\deltadist\R_deltadist_1_5_gldrs.txt'
# from papers/dhk/2019/figures
# copied from sep4_deltadist/deltadist
# oct10/deltadist
# file was called 'R_variances_1_4_gldrs.txt'

library(fishmethods)

Fig9 <- function(
AMLR.area = c('SA'),
n.rep=9,
NASC.yrs = c(2001:2009,2011),
n.gldrs = c(1,4,5)) {

  shiptable <- as.matrix(read.table(
              paste('tables/shiptable',AMLR.area,'.txt',sep=''),header=TRUE),
              header = TRUE)
  ship.s <- shiptable[,'sum'] # NASC, need to convert to biom
  ship.sds <- shiptable[,'sd']

  gldr.1.list <- gldr.4.list <- gldr.5.list <- list()
  for (iyr in 1:length(NASC.yrs)){
    gldr.1.list[[iyr]] <- as.matrix(read.table(paste('tables/yo_sums/',n.gldrs[1],
                   '_',AMLR.area,'_',NASC.yrs[iyr],'_150_yo_sums.txt',sep=''),header=TRUE))
    gldr.4.list[[iyr]] <- as.matrix(read.table(
                   paste('tables/yo_sums/',n.gldrs[2],'_',AMLR.area,'_',NASC.yrs[iyr],
                   '_150_yo_sums.txt',sep=''),header=TRUE))
    gldr.5.list[[iyr]] <- as.matrix(read.table(
                   paste('tables/yo_sums/',n.gldrs[3],'_',AMLR.area,'_',NASC.yrs[iyr],
                   '_150_yo_sums.txt',sep=''),header=TRUE))
    }

  delta.dist.1 <- delta.dist.4 <-  delta.dist.5 <- mean.sd.1 <- mean.sd.4 <-
    array(dim=c(length(NASC.yrs),3,n.rep)) # replicates of a NASC.yrs X 3 stats array
  dimnames(delta.dist.1) <- dimnames(delta.dist.4) <- dimnames(delta.dist.5) <- 
    list(NASC.yrs,c('ship','mean','var'))
  dimnames(mean.sd.1) <- dimnames(mean.sd.4) <-
    list(NASC.yrs,c('ship','mean','sd'))
  for (iyr in 1:length(NASC.yrs))
    for(i.rep in 1: n.rep){
      delta.dist.1[iyr,,i.rep] <- c(ship.s[iyr],deltadist(gldr.1.list[[iyr]][,i.rep]))
      delta.dist.4[iyr,,i.rep] <- c(ship.s[iyr],deltadist(gldr.4.list[[iyr]][,i.rep]))
      delta.dist.5[iyr,,i.rep] <- c(ship.s[iyr],deltadist(gldr.5.list[[iyr]][,i.rep]))
      mean.sd.1[iyr,,i.rep] <- c(ship.s[iyr],c(mean(gldr.1.list[[iyr]][,i.rep]),
                       sd(gldr.1.list[[iyr]][,i.rep])
                      ))
      mean.sd.4[iyr,,i.rep] <- c(ship.s[iyr],c(mean(gldr.4.list[[iyr]][,i.rep]),
                       sd(gldr.4.list[[iyr]][,i.rep])
                      ))
    }
log.m <-log(apply(delta.dist.1[,'mean',],1,mean,na.rm=TRUE))
log.sd <-log(apply(delta.dist.1[,'mean',],1,sd,na.rm=TRUE))

#plt.name <- paste('delta_dist_1gldr_',n.rep,'_reps.pdf',sep='')
#pdf(file = plt.name) #,width=24,height=18)
par(mfrow=c(1,1),cex=1.5,cex.lab=1.5,oma=c(2,2,2,0))
y.lim <- c(min(log(ship.s),log.m/
               (log.sd/log.m)),
           max(log(ship.s),log.m/
               (log.sd/log.m))
          )
plot(NASC.yrs,log(ship.s),cex=1.5,pch=19,lwd=3,ylim=y.lim,
    main = paste('Delta.dist 1 glider_',AMLR.area,'_',n.rep,' reps',sep=''),
           ylab='',xlab='Year')
title(,ylab=expression(paste('Mean log(',s[a],')')),mgp=c(2.5,1,0))
lines(NASC.yrs,log(apply(delta.dist.1[,'mean',],1,mean,na.rm=TRUE)),lwd=3)
segments(
         NASC.yrs,log.m/(log.sd/log.m),
         NASC.yrs,log.m*(log.sd/log.m),
	 lwd=3
        ) 
#dev.off()
} # end Fig9()
