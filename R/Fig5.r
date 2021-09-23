# modified from 'C:\zot\glider\2021\9sep\sep7_mac\rmarkdown'

Fig5 <- function(NASC.yrs = c(2001:2009,2011),
  AMLR.area = "SA",
  leg = 1,n.rep = 9,n.gldr =c(1),
  qntl.vals = c(0.97,0.98,0.99,0.999,1,"sum","sd"),
  depths = c(150,1000),azfp.off = c(150,150)
  ) {
  gldr.smpls.m <- gldr.smpls.sd <- vector()
  for(ipth in 1:length(n.gldr))
    for(idpth in 1:length(depths)){
    yo.count <- read.table(paste('tables/',n.gldr[ipth],'_', AMLR.area,
                '_',depths[idpth],"_yo_count.txt",sep=''))
    ship.sums <- as.matrix(read.table('tables/shiptable.txt',header=TRUE)[,2:8])
    gldr.smpls.yrs <- array(dim=c(n.rep,length(qntl.vals),length(NASC.yrs)))
    dimnames(gldr.smpls.yrs) <- list(1:n.rep,c(qntl.vals),NASC.yrs)
    mae.table <- rmse.table <-array(dim=c(length(NASC.yrs),length(qntl.vals)))
    dimnames(mae.table) <- dimnames(rmse.table) <- list(NASC.yrs,c(qntl.vals))
    for(iyr in 1:length(NASC.yrs)){
      gldr.smpls <- read.table(paste('tables/gldrtables/',n.gldr[ipth],'_', AMLR.area,'_',NASC.yrs[iyr],
                  "_",depths[idpth],"m_gldrtable.txt",sep=""),header=TRUE)
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
          rmse.tmp[i,] <- as.numeric(sqrt(((unlist(gldr.smpls[i,,iyr])-ship.sums[iyr,])^2)))
          }
        }
      gldr.smpls.m[iyr] <- mean(gldr.smpls.yrs[,"sum",iyr])
      gldr.smpls.sd[iyr] <- sd(gldr.smpls.yrs[,"sum",iyr])
      mae.table[iyr,] <- apply(mae.tmp,2,sum) 
      rmse.table[iyr,] <- sqrt(apply(rmse.tmp,2,sum)/
                                 nrow(gldr.smpls.yrs[,,iyr])) 
    } # end iyr
    mae.table <- cbind(mae.table, yo.count = as.numeric(yo.count))
    rmse.table <- cbind(rmse.table, yo.count = as.numeric(yo.count))
    dimnames(yo.count) <- list("yo.count",NASC.yrs)
    write.table(round(mae.table,1),paste(n.gldr[ipth],"_annual_mae.txt",sep=''))
    write.table(round(rmse.table,1),paste(n.gldr[ipth],"_annual_rmse.txt",sep=''))
 #  plt.name <- "Glider & Ship Sv_log.pdf"
 #   pdf(file = plt.name) #,width=24,height=18)
    par(mfrow=c(1,1),cex=1.3,cex.lab=1.3,oma=c(2,2,2,0))
    y.lim <- c(min(c(log(ship.sums[,"sum"]),
                   log(gldr.smpls.yrs[,"sum",][gldr.smpls.yrs[,"sum",]>0]),na.rm=TRUE)),
             max(c(log(ship.sums[,"sum"]),
                   log(gldr.smpls.yrs[,6,])),na.rm=TRUE))
    plot(NASC.yrs,log(ship.sums[,"sum"]),type="l",col="red",
      ylim=y.lim,lwd=3,main=paste(depths[idpth],"m max yo depth",sep=""),
          ylab="",xlab = "Year")
      title(ylab = expression(paste("log(s" [a],")")), mgp=c(2.5,1,0))
    for(i.rep in 1:n.rep){
      lines(NASC.yrs,log(gldr.smpls.yrs[i.rep,"sum",]),col="gray",lwd=3)
      points(NASC.yrs,log(ship.sums[,"sum"]),col="red",pch=19)
      }
    lines(NASC.yrs,log(gldr.smpls.m),col="blue",lwd=3)
    log.cv <- exp(sqrt(log(1+(gldr.smpls.sd/gldr.smpls.m)^2)))
    segments(NASC.yrs,log(gldr.smpls.m * log.cv),
           NASC.yrs,log(gldr.smpls.m / log.cv),
           col="blue",lwd=3)
    points(NASC.yrs,log(gldr.smpls.m),
         pch=19,col="blue")
    lines(NASC.yrs,log(ship.sums[,"sum"]),col="red",lwd=3)
  } # end ipth, idpth
} # end function(Fig5)







