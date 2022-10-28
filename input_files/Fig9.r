# calculate delta-distribution estimates of mean and SD
# of glider samples for comparison to logarithmic CV
# method used in this glider simulation

Fig9 <- function(NASC.yrs = c(2001:2009,2011),
  AMLR.areas = c('SA','WA'),
  depths = c(
            150 #,200,300,400,500,700,1000
           ), 
  n.gldr = c('1') # one set of glider combinations at a time
  ){

  if(!dir.exists('Fig9'))
    dir.create('Fig9')
  for(i.gldr in 1:length(n.gldr))
  for(i.depth in 1:length(depths)){  
  #for(i.gldr in c(1,3,5)) # 1, 3 and 5 gliders only
  #for(i.depth in c(1,7)){ # 150 m and 1000 m only)  
  #F8.depth <- depths[i.depth]
  for(i.area in 1: length(AMLR.areas)){
    #path.F8 <- paste("glider/2021/1jan/jan12/",AMLR.areas[i.area],"/",
    #                 n.gldr[i.gldr],'gldrs/',sep="")
    ship.stats <- read.table(paste('tables/shiptable',AMLR.areas[i.area],'.txt',sep=''),header=TRUE)[,2:8]
    ship.m <- ship.stats[,'sum']
    ship.sd <- ship.stats[,'sd']

      gldr.m.l <- gldr.sd <- gldr.m  <- list() # dim for each year will be n.yos X n.reps
      gldr.t.m <- gldr.t.sd <- vector() # the mean or sd of all yos for each year
      for(iyr in 1:length(NASC.yrs)){
        gldr.m.l[[iyr]] <- read.table(paste('tables/yo_sums/',n.gldr[i.gldr],'_',AMLR.areas[i.area],'_', 
               NASC.yrs[iyr],'_',depths[i.depth],'_yo_sums.txt',sep=''),header=TRUE)
               # only works for one cglider at a time
        gldr.t.m[iyr] <- mean(unlist(gldr.m.l[[iyr]]),na.rm=TRUE)
        gldr.t.sd[iyr] <- sd(unlist(gldr.m.l[[iyr]]),na.rm=TRUE)
        }
      n.reps <- ncol(gldr.m.l[[iyr]]) # yos (given depth) X replicates

      gldr.m <- gldr.sd <-  array(dim=c(length(NASC.yrs),n.reps))
      for(iyr in 1:length(NASC.yrs))
        for (i.rep in 1:n.reps){
        gldr.m[iyr,i.rep] <- mean(gldr.m.l[[iyr]][,i.rep],na.rm=TRUE)
        gldr.sd[iyr,i.rep] <- sd(gldr.m.l[[iyr]][,i.rep],na.rm=TRUE)
        }
    names(gldr.m) <- names(gldr.sd) <- NASC.yrs

    # test different glider years for coverage of 68% for the SD.
    coverage.p <- vector()
    for(iyr in 1:length(NASC.yrs)){
      gldr.cv <- exp(sqrt(log(1+(apply(gldr.m,2,sd,na.rm=TRUE)/
                            apply(gldr.m,2,mean,na.rm=TRUE))^2)))
      coverage.p[iyr] <- table(ship.m[iyr] < gldr.m[iyr,]*gldr.cv & 
        ship.m[iyr] > gldr.m[iyr,]/gldr.cv)['TRUE']/
        length(gldr.m[iyr,]) # proportion inside nominal SE
      }
    names(coverage.p) <- NASC.yrs
    write.table(coverage.p,paste(getwd(),'/Fig9/','coverage_p_',
             n.gldr[i.gldr],'_',AMLR.areas[i.area],'_',depths[i.depth],'m.txt',sep=''))
  } # end i.area
  combo.files <-c(paste('coverage_p_',n.gldr[i.gldr],'_SA_',
      depths[i.depth],'m.txt',sep=''),
      paste('coverage_p_',n.gldr[i.gldr],'_WA_',depths[i.depth],'m.txt',sep=''))
  combo.table <- read.table(paste('Fig9/',combo.files[1],sep=''))
  for (i in 2: length(combo.files))
    combo.table <- cbind(combo.table,read.table(paste('Fig9/',combo.files[i],sep='')))
  colnames(combo.table) <- c(paste('BS_',depths[i.depth],'m',sep=''),
                           paste('CS_',depths[i.depth],'m',sep=''))
  combo.table <- rbind(combo.table[1:9,],c(NA,NA),combo.table[10,])
  rownames(combo.table) <- 2001:2011

  #plt.name <- paste('Fig9/Fig9_',n.gldr[i.gldr],'gldr_',depths[i.depth],'m_reps.jpeg',sep='')
  #jpeg(file = plt.name) #,width=24,height=18)
  par(cex=1.3)
  #y.lim <- c(min(combo.table,0.68,na.rm=TRUE),max(combo.table,na.rm=TRUE))
  y.lim <- c(0,1)
  plot(rownames(combo.table),combo.table[,1],type='p',pch=19,lty=1,lwd=3,col='blue',
     main=paste('Coverage probabilities',substr(n.gldr[i.gldr],1,5),depths[i.depth],'m'),
     ylim=y.lim,ylab='Coverage',xlab='Year')

  points(rownames(combo.table),combo.table[,2],pch=19,lty=2,lwd=3,col='red')
  lines(rownames(combo.table),combo.table[,1],lty=1,lwd=3,col='blue')
  lines(rownames(combo.table),combo.table[,2],lty=2,lwd=3,col='red')
  legend(2008,0.4,c(paste('BS ',sep=''),
                 paste('CS ',sep='')),
       lty=c(1,2),lwd=3,col=c('blue','red'),cex=0.85)
  abline(h=0.68, lty = 2, lwd = 2)
  #graphics.off()
  write.table(round(combo.table,2),paste('Fig9/yo_coverage_table_',depths[i.depth],'m.txt',sep=''))
  } # end i.gldr, i.dpth

} # end Fig9()
