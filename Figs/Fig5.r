Fig5 <- function(NASC.yrs = c(2001:2009,2011),
  AMLR.area = "SA",
  leg = 1,n.rep = 9,
  qntl.vals = c(0.97,0.98,0.99,0.999,1,"sum","sd"),
  depths = c(150)
  ) {
  gldr.smpls.m <- gldr.smpls.sd <- vector()
#  for(i.depth in 1:length(depths)){
    yo.count <- read.table("yo_count.txt")
    ship.sums <- as.matrix(read.table("shiptable.txt",header=TRUE)[,2:8])
    if(leg==1)
      UFFS <- read.csv("UFFS_leg1.csv",sep=",")
    if(leg==2)
      UFFS <- read.csv("UFFS_leg2.csv",sep=",")
    UFFS.yrs <- NASC.yrs - (1995+1)
    UFFS.dat <- UFFS[UFFS.yrs,AMLR.area] # find UFFS yrs and area
    gldr.smpls.biom <- gldr.smpls.yrs <- array(dim=c(n.rep,length(qntl.vals),length(NASC.yrs)))
    dimnames(gldr.smpls.biom) <- dimnames(gldr.smpls.yrs) <- list(1:n.rep,c(qntl.vals),NASC.yrs)
    mae.table <- rmse.table <-array(dim=c(length(NASC.yrs),length(qntl.vals)))
    dimnames(mae.table) <- dimnames(rmse.table) <- list(NASC.yrs,c(qntl.vals))
    ship.sums.biom <- array(dim=c(length(NASC.yrs),7))
    dimnames(ship.sums.biom) <- list(NASC.yrs,c(qntl.vals))
    for(iyr in 1:length(NASC.yrs)){
      gldr.smpls <- read.table(paste('gldrtables/',NASC.yrs[iyr],"_",depths,
                  "m_gldrtable.txt",sep=""),header=TRUE)
      gldr.smpls.biom[,,iyr] <- as.matrix(UFFS.dat[iyr] * gldr.smpls)
      gldr.smpls.yrs[,,iyr] <- as.matrix(gldr.smpls)
      ship.sums.biom[iyr,] <- as.numeric(UFFS.dat[iyr]*ship.sums[iyr,]) 
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
    write.table(round(mae.table,1),"annual_mae.txt")
    write.table(round(rmse.table,1),"annual_rmse.txt")
} # end function(Fig5)







