
library(fishmethods)

Fig9 <- function(
  AMLR.area = c('SA'),
  n.rep=9,
  NASC.yrs = c(2001:2009,2011),
  n.gldr = c(1,4,5)) {

  for(i.area in 1:length(AMLR.area)){
    shiptable <- as.matrix(read.table(
              paste('tables/shiptable',AMLR.area[i.area],'.txt',sep=''),header=TRUE),
              header = TRUE)
    ship.s <- shiptable[,'sum'] # NASC, need to convert to biom
    ship.sds <- shiptable[,'sd']

    gldr.list <- delta.dist <- mean.sd <- list()
    for (i.gldr in 1:length(n.gldr)){
      gldr.list[[i.gldr]] <- list()
      for (iyr in 1:length(NASC.yrs)){
        gldr.list[[i.gldr]][[iyr]] <- as.matrix(read.table(paste('tables/yo_sums/',n.gldr[i.gldr],
                   '_',AMLR.area[i.area],'_',NASC.yrs[iyr],'_150_yo_sums.txt',sep=''),header=TRUE))
        } # end iyr

      delta.dist[[i.gldr]] <- mean.sd[[i.gldr]] <-
        array(dim=c(length(NASC.yrs),3,n.rep)) # replicates of a NASC.yrs X 3 stats array
      dimnames(delta.dist[[i.gldr]]) <- 
        list(NASC.yrs,c('ship','mean','var'))
      dimnames(mean.sd[[i.gldr]]) <- 
        list(NASC.yrs,c('ship','mean','sd'))
      for (iyr in 1:length(NASC.yrs))
        for(i.rep in 1: n.rep){
          delta.dist[[i.gldr]][iyr,,i.rep] <- c(ship.s[iyr],deltadist(gldr.list[[i.gldr]][[iyr]][,i.rep]))
          mean.sd[[i.gldr]][iyr,,i.rep] <- c(ship.s[iyr],c(mean(gldr.list[[i.gldr]][[iyr]][,i.rep]),
                       sd(gldr.list[[i.gldr]][[iyr]][,i.rep])
                      ))

        } # end iyr, i.rep
      } # end i.gldr

  log.m <-log(apply(delta.dist[[1]][,'mean',],1,mean,na.rm=TRUE))
  log.sd <-log(apply(delta.dist[[1]][,'mean',],1,sd,na.rm=TRUE))

  #plt.name <- paste('delta_dist_1gldr_',n.rep,'_reps.pdf',sep='')
  #pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,cex.lab=1.5,oma=c(2,2,2,0))
  y.lim <- c(min(log(ship.s),log.m/
               (log.sd/log.m)),
           max(log(ship.s),log.m/
               (log.sd/log.m))
          )
  plot(NASC.yrs,log(ship.s),cex=1.5,pch=19,lwd=3,ylim=y.lim,
    main = paste('Delta.dist 1 glider_',AMLR.area[i.area],'_',n.rep,' reps',sep=''),
           ylab='',xlab='Year')
  title(,ylab=expression(paste('Mean log(',s[A],')')),mgp=c(2.5,1,0))
  lines(NASC.yrs,log(apply(delta.dist[[1]][,'mean',],1,mean,na.rm=TRUE)),lwd=3)
  segments(
         NASC.yrs,log.m/(log.sd/log.m),
         NASC.yrs,log.m*(log.sd/log.m),
	 lwd=3
        ) 
  #dev.off()
  } # end i.area
} # end Fig9()
