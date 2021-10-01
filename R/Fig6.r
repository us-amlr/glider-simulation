# uses large cgldr matrix of combined replicates
# from 'aug20/Fig6/orig_contour_code/R_filled_contours.txt'
# combined with 'aug25/R_read_NAS_yrs.txt' and
# 'aug25/R_5_contour.txt'
# got rid of 'gldr.depth' cbind calculations
# replaced 'gldrXyrXdepth' with 'gldr.smpls'

strt <- Sys.time()


Fig6 <- function(depths = c(150), AMLR.area = 'SA', max.z.b = 50, NASC.yrs = c(2001:2009,2011), n.reps = 9, n.gldr = c(1) ){
  dat.biom <- list() # dat.NASC also a list in 'R_5_contour_cgldr_mac.txt'
  cgldr.pth <- paste(getwd(),'/tables/cgldrs/',sep='')
  s.pth <- paste(getwd(),'/tables/smpl_paths/',sep='')
  dbins <- rev(c(1:max.z.b)) # 5m depth bins sampled by gliders
  depth.means <- array(dim=c(length(NASC.yrs),max.z.b))
  ship.NASC <- list() # dat.NASC[[iyr]][[i.depth]] in 'R_5_contour_cgldr_mac.txt'
  for (iyr in 1:length(NASC.yrs)){
    ship.NASC[[iyr]] <- read.csv(paste('tables/NASC_yrs/',
      n.gldr,'_',AMLR.area,'_',NASC.yrs[iyr],'_NASC_DAT.csv',sep=''))
    ship.NASC[[iyr]] <- ship.NASC[[iyr]][,-1] # drop first column of
                                              # Process_ID and Interval values
    depth.means[iyr,] <- apply(ship.NASC[[iyr]],2,mean,na.rm=TRUE)
    }

  #######################

  # Glider samples
  # read cgldr tables and combine the replicates
  # calculate annual and max yo depth glider NASC in depth bins
  # (from 'R_5_contour_cgldr_mac.txt')

  gldr.smpls <- gldr.smpls.1 <-gldr.smpls.5 <- array(dim=c(max.z.b,length(NASC.yrs),n.gldr))
  c.gldrs  <- list()
  for(i.gldr in 1:n.gldr){
    c.gldrs[[i.gldr]] <- list()
    for(iyr in 1:length(NASC.yrs)){
     # 'depths' is max yo depth in meters
     tmp <- read.table(paste(cgldr.pth,n.gldr,'_', AMLR.area,'_',NASC.yrs[iyr]
               ,'_',depths,'_cgldrs.txt',sep=''),header=TRUE)
               # only works for one cglider so far
        c.gldrs[[i.gldr]][[iyr]] <- array(dim=c(dim(tmp)))
        c.gldrs[[i.gldr]][[iyr]] <- tmp
        #cols.rep[iyr] <- ncol(c.gldrs[[i.gldr]][[iyr]])/n.reps
        smpl.pths <- read.table(paste(s.pth,n.gldr,'_', AMLR.area,'_',NASC.yrs[iyr],'_',
               depths,'_spths.txt',sep='')) # smpl.pths is for 1 replicate
        mask <- array(NA,dim=c(nrow(smpl.pths),
              (ncol(smpl.pths)*n.reps)))
	mask <- as.data.frame(mask)
	rep.cols <- c(0,seq(1:n.reps) * ncol(smpl.pths))
	for (i.rc in 1:(length(rep.cols)-1)){
	  mask[,(rep.cols[i.rc]+1):rep.cols[i.rc+1]] <- smpl.pths # dimension of mask turns NULL
	  }

        tmp2 <- c.gldrs[[i.gldr]][[iyr]][mask>0]
	dim(tmp2) <- dim(c.gldrs[[i.gldr]][[iyr]])
        gldr.smpls[,iyr,i.gldr] <- apply(tmp2,1,mean,na.rm=TRUE)
	gldr.smpls.1[,iyr,i.gldr] <- apply(tmp2[,1:rep.cols[2]],1,mean,na.rm=TRUE)
	gldr.smpls.5[,iyr,i.gldr] <- apply(tmp[,1:rep.cols[6]],1,mean,na.rm=TRUE)
    } # end iyr
  } # end i.gldr
	

  #######################
  #plt.name <- 'contour_ship.pdf'
  #pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,cex.lab=1.5,oma=c(2,2,2,0))
  depth.means[is.na(depth.means)]=0
  # color plot from 'aug20/fig6/orig_code/R_filled_contours.txt'
  filled.contour(x=c(1:length(NASC.yrs)),y=c(1:length(dbins)),
        z=depth.means[,50:1],
        axes = FALSE,#key.axis=FALSE,
        main=paste('Population densities'),cex.main=1.5,
        plot.axes = {
          axis(1,at = 1:length(NASC.yrs),labels=NASC.yrs)
	  axis(2,at = length(dbins):1,labels=rev(dbins)*5)
          },
        key.axes = log(axis(4)),
	key.title = title(expression(paste('(',s[a],')'))),
        col=terrain.colors(11), #viridis(10),
        xlab='Year',ylab='Depth (m)',
        levels=c(0,0.1,0.5,1,2,3,5,8,12,24,50,max(depth.means,na.rm=TRUE))
        )
#  dev.off()

  #######################
#  plt.name <- paste(path.write,dbins,'contour_',n.reps,'reps.pdf',sep='')
#  pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,cex.lab=1.5,oma=c(2,2,2,0))
  gldr.smpls[is.na(gldr.smpls)]=0
  filled.contour(x=c(1:length(NASC.yrs)),y=c(1:length(dbins)),
        z=t(gldr.smpls[50:1,1:ncol(gldr.smpls),i.gldr]),
        axes = FALSE,#key.axis=FALSE,
        main=paste('Mean densities from\n' ,n.reps,' replicate gliders',sep=''),cex.main=1.5,
        plot.axes = {
          axis(1,at = 1:length(NASC.yrs),labels=NASC.yrs)
	  axis(2,at = length(dbins):1,labels=rev(dbins)*5)
          },
        key.axes = log(axis(4)),
	key.title = title(expression(paste('(',s[a],')'))),
        col=terrain.colors(11), #viridis(10),
        xlab='Year',ylab='Depth (m)',
        levels=c(0,0.1,0.5,1,2,3,5,8,12,24,50,
	       max(gldr.smpls,na.rm=TRUE))
        )
#  dev.off()

  #######################
#  plt.name <- paste(path.write,dbins,'contour_1rep.pdf',sep='')
#  pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,cex.lab=1.5,oma=c(2,2,2,0))
  gldr.smpls[is.na(gldr.smpls)]=0
  filled.contour(x=c(1:length(NASC.yrs)),y=c(1:length(dbins)),
        z=t(gldr.smpls.1[50:1,1:ncol(gldr.smpls.1),i.gldr]),
        axes = FALSE,#key.axis=FALSE,
        main=paste('Mean densities from\n' ,'1 glider',sep=''),cex.main=1.5,
        plot.axes = {
          axis(1,at = 1:length(NASC.yrs),labels=NASC.yrs)
	  axis(2,at = length(dbins):1,labels=rev(dbins)*5)
          },
        key.axes = log(axis(4)),
	key.title = title(expression(paste('(',s[a],')'))),
        col=terrain.colors(11), #viridis(10),
        xlab='Year',ylab='Depth (m)',
        levels=c(0,0.1,0.5,1,2,3,5,8,12,24,50,
	       max(gldr.smpls.1,na.rm=TRUE))
        )
#  dev.off()

  #######################
#  plt.name <- paste(path.write,dbins,'contour_5reps.pdf',sep='')
#  pdf(file = plt.name) #,width=24,height=18)
  par(mfrow=c(1,1),cex=1.5,cex.lab=1.5,oma=c(2,2,2,0))
  gldr.smpls[is.na(gldr.smpls)]=0
  filled.contour(x=c(1:length(NASC.yrs)),y=c(1:length(dbins)),
        z=t(gldr.smpls.5[50:1,1:ncol(gldr.smpls.5),i.gldr]),
        axes = FALSE,#key.axis=FALSE,
        main=paste('Mean densities from\n' ,'5_replicate gliders',sep=''),cex.main=1.5,
        plot.axes = {
          axis(1,at = 1:length(NASC.yrs),labels=NASC.yrs)
	  axis(2,at = length(dbins):1,labels=rev(dbins)*5)
          },
        key.axes = log(axis(4)),
	key.title = title(expression(paste('(',s[a],')'))),
        col=terrain.colors(11), #viridis(10),
        xlab='Year',ylab='Depth (m)',
        levels=c(0,0.1,0.5,1,2,3,5,8,12,24,50,
	       max(gldr.smpls.5,na.rm=TRUE))
        )
#  dev.off()

  strt
  Sys.time()
} # end function
