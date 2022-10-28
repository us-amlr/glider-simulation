library(viridis)

strt <- Sys.time()

Fig6 <- function(depths = c(150), AMLR.area = 'SA', max.z.b = 50,
                 NASC.yrs = c(2001:2009,2011), n.rep = 9, n.gldr = n.gldr ){
  dat.biom <- list() 
  cgldr.pth <- paste(getwd(),'/tables/cgldrs/',sep='')
  s.pth <- paste(getwd(),'/tables/smpl_paths/',sep='')
  dbins <- rev(c(1:max.z.b)) # 5m depth bins sampled by gliders
  depth.means <- array(dim=c(length(NASC.yrs),max.z.b))
  ship.NASC <- list() 
  for (iyr in 1:length(NASC.yrs)){
    ship.NASC[[iyr]] <- read.csv(paste('tables/NASC_yrs/',
      n.gldr[1],'_',AMLR.area,'_',NASC.yrs[iyr],'_NASC_DAT.csv',sep=''))
    ship.NASC[[iyr]] <- ship.NASC[[iyr]][,-1] # drop first column of
                                              # Process_ID and Interval values
    depth.means[iyr,] <- apply(ship.NASC[[iyr]],2,mean,na.rm=TRUE)
    }
  ####################### 
  # Plot the population values the gliders are sampling
  # by strata, depth bin and year
  #plt.name <- 'contour_ship.pdf'
  #pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=0.9,cex.lab=1.5,oma=c(2,2,2,0))
  depth.means[is.na(depth.means)]=0
  # color plot from 'aug20/fig6/orig_code/R_filled_contours.txt'
  filled.contour(x=c(1:length(NASC.yrs)),y=c(1:length(dbins)),
        z=depth.means[,50:1],
        axes = FALSE,#key.axis=FALSE,
        main=paste('Population densities',sep=''),cex.main=1.5,
        plot.axes = {
          axis(1,at = 1:length(NASC.yrs),labels=NASC.yrs)
	  #axis(2,at = length(dbins),labels=rev(dbins)*5)
	  axis(2,at = seq(1,max(dbins),by=4),
	         labels=(rev(seq(1,max(dbins),by=4)*5)+5))
          contour(x=c(1:length(NASC.yrs)),y=c(1:length(dbins)),
                 z=depth.means[,50:1],add=TRUE,labcex=2,
                 levels=c(0,0.5,2,5,12,24,50,75,100,150,
                 max(depth.means,na.rm=TRUE)),col='white')
          },
        key.axes = log(axis(4)),
	key.title = title(expression(paste('(',s[A],')'))),
        color.palette =  viridis, # col = terrain.colors(11),
        xlab='Year',ylab='Depth (m)',
        levels=c(0,0.5,2,5,12,24,50,75,100,150,
               max(depth.means,na.rm=TRUE))
        )
#  dev.off()

  #######################
  # Glider samples:
  # read the cgldr data (the regions being sampled for all replicate dives of
  # 1-5 combined gliders sampling from the surface down to 50 5m depth bins) 
  # and combine the replicates to calculate annual NASC sampled by the gliders.
  # 
  # 'tmp' is the matrix being sampled by each replicate. 
  # dim(tmp) = (max.z.b, ncol(smpl.pths * n.rep)).
  # 'mask' locates the glider samples in each dive through the cgldr regions.
  # 'tmp2' is the matrix of samples obtained by the glider(s) for each replicate.
  # 'rep.cols' are the starting columns of 'tmp2' for each replicate.

  gldr.reps <- array(dim=c(max.z.b,length(NASC.yrs),length(n.gldr)))
  c.gldrs  <- list()
  for(i.gldr in 1:length(n.gldr)){
    c.gldrs[[i.gldr]] <- list()
    for(iyr in 1:length(NASC.yrs)){
     # 'depths' is max yo depth in meters
     tmp <- read.table(paste(cgldr.pth,n.gldr[i.gldr],'_', AMLR.area,'_',NASC.yrs[iyr]
               ,'_',depths,'_cgldrs.txt',sep=''),header=TRUE)
               # only works for one cglider at a time
        c.gldrs[[i.gldr]][[iyr]] <- array(dim=c(dim(tmp)))
        c.gldrs[[i.gldr]][[iyr]] <- tmp
        smpl.pths <- read.table(paste(s.pth,n.gldr[i.gldr],'_', AMLR.area,'_',NASC.yrs[iyr],'_',
               depths,'_spths.txt',sep='')) # smpl.pths is for 1 replicate
        mask <- array(NA,dim=c(nrow(smpl.pths),
              (ncol(smpl.pths)*n.rep)))
	mask <- as.data.frame(mask)
	rep.cols <- c(0,seq(1:n.rep) * ncol(smpl.pths))
	for (i.rc in 1:(length(rep.cols)-1)){
	  mask[,(rep.cols[i.rc]+1):rep.cols[i.rc+1]] <- smpl.pths 
	  }
        tmp2 <- c.gldrs[[i.gldr]][[iyr]][mask>0]
	dim(tmp2) <- dim(mask) 
        gldr.reps[,iyr,i.gldr] <- apply(tmp2,1,mean,na.rm=TRUE) # mean of all replicates
    } # end iyr
	

  #######################
#  plt.name <- paste('gldr_contour_',n.rep,'reps.pdf',sep='')
#  pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=0.9,cex.lab=1.5,oma=c(2,2,2,0))
  gldr.reps[is.na(gldr.reps)]=0
  #dev.new()
  filled.contour(x=c(1:length(NASC.yrs)),y=c(1:length(dbins)),
        z=as.matrix(t(gldr.reps[50:1,,i.gldr])),
        axes = FALSE,#key.axis=FALSE,
        main=ifelse(n.gldr[i.gldr]==1,paste('Mean densities from\n' ,n.gldr[i.gldr],' glider',sep=''),
                    paste('Mean densities from\n' ,n.gldr[i.gldr],' gliders',sep='')),
        cex.main=1.5,
        plot.axes = {
          axis(1,at = 1:length(NASC.yrs),labels=NASC.yrs)
	  #axis(2,at = length(dbins):1,labels=rev(dbins)*5)
	  axis(2,at = seq(1,max(dbins),by=4),
	         labels=(rev(seq(1,max(dbins),by=4)*5)+5))
          contour(x=c(1:length(NASC.yrs)),y=c(1:length(dbins)),
                 z=as.matrix(t(gldr.reps[50:1,,i.gldr])),add=TRUE,labcex=2,
                 levels=c(0,0.5,2,5,12,24,50,75,100,150,
                       max(depth.means,na.rm=TRUE)),col='white')
          },
        key.axes = log(axis(4)),
	key.title = title(expression(paste('(',s[A],')'))),
        color.palette =  viridis, # col = terrain.colors(11),
        xlab='Year',ylab='Depth (m)',
        levels=c(0,0.5,2,5,12,24,50,75,100,150,
                       max(depth.means,na.rm=TRUE))
        )

  #  dev.off()
  } # end i.gldr
} # end function