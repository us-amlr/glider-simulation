# autonomous underwater glider sampling simulation
# calls gldry()

gldrs <- function(NASC.yrs = c(2001:2009,2011),AMLR.area = 'SA',n.rep = 9,n.gldr = c(1,2,3,4,5),
                  save.tables = 1,max.NASC.m = 250,depths = c(150,200,300,400,500,700,1000),azfp.off = 
                  c(150,150,150,150,150,150,150),
                  qntl.vals = c(0.97,0.98,0.99,0.999,1),smpl.st = 0 
                  ){
max.z.m <- max.NASC.m #??
if(save.tables == 1)
  if(!dir.exists('tables')){
    dir.create('tables')
    dir.create('tables/yo_sums')
    dir.create('tables/gldrtables')
    dir.create('tables/gldr_strt')
    dir.create('tables/NASC_yrs') # operating model bin densities
    dir.create('tables/cgldrs') # glider sampled bin densities
    dir.create('tables/smpl_paths/') # bins sampled by gliders
    dir.create('tables/gldr_stats/') # annual log means and SDs of glider samples
    dir.create('tables/yo_count/') # number of glider yos/year
    }

ship.Sv <- vector() # ship acoustic mean integrated returns each year
# load the databases
#load('NASC_leg1.RData') # acoustic data to be sampled, 199772 100m long X 100 5m deep cells
#leg1 <- read.csv('leg1.csv') # date and sampling strata for each of 199772 columns of acoustic data
for(i.gldr in n.gldr) # i.gldr indexes the number of gliders that are sampling at the same time
for(i.depth in 1:length(depths)){
  std.yo.depth.m <- depths[i.depth] # max yo depth
  shutoff.azfp.m <- azfp.off[i.depth]

  # assign the area, year, and depths to be sampled
  yo.count <- vector()
  ship.table <- array(dim=c(length(NASC.yrs),length(qntl.vals)+3))
  colnames(ship.table) <- c('year',c(qntl.vals,'sum','sd'))
  options(stringsAsFactors=FALSE)

  gldr.sums <- gldr.sds <- array(dim=c(n.rep,length(NASC.yrs)))
  dimnames(gldr.sums) <- dimnames(gldr.sds) <- list(paste('r.',c(1:n.rep),sep=''),NASC.yrs)
  ship.sums <- vector(length=length(NASC.yrs))
  names(ship.sums)  <- NASC.yrs
  for(iyr in 1:length(NASC.yrs)){
    gldr.strt.table <- gldr.strt.table.nm <- array(dim=c(i.gldr,n.rep))
    smpl.yr <- c(NASC.yrs[iyr])
    NASC.vals <- leg1[substr(leg1$Date_M,1,4) == smpl.yr &
                    #leg1$Depth_mean <= shutoff.azfp.m+100 &  # -dhk
		    leg1$area == AMLR.area,]
    # NASC.vals is the subset of the entire ship database
    # that will be sampled with the glider. It is used to
    # find the Process_ID and Interval values for sampling
    #NASC.vals <- NASC.vals[is.na(NASC.vals)] <- 0 #-dhk

    yo.locs <- PID <- vector()  
    # find the process and interval IDs (column names) in NASC.vals
    PID.uniq <- unique(NASC.vals$Process_ID)
    k <- 0
    for (i in 1:length(PID.uniq)){
      PID[i] <- PID.uniq[i]
      PID.int <- unique(NASC.vals$Interval[NASC.vals$Process_ID == PID[i]])
      for(j in 1:length(PID.int)){
        k <- k+1
        yo.locs[k] <- paste(PID.uniq[i],'.',PID.int[j],sep='')
        }
      # yo.locs are the column names to be sampled in NASC.leg1
      # NASC.leg1 is the ship NASC matrix from shipNASC.RData
      }
    # use 'gldry()' outputs (out$z.bin, out$d.bin, out$dive.id)
    # to calculate the glider locations (smpl.pth) in the NASC grid
    #source(paste(getwd(),'/glider_yos.r',sep=''))
    source('gldry.r')
    out.list <- gldry(std.yo.depth.m,shutoff.azfp.m,max.z.m,yo.locs)
    n.yos <- out.list[[1]]
    out <- out.list[[2]]
    azfp.range.m <- out.list[[3]]
    bin.depth.m <- out.list[[4]]
    yo.count[iyr] <- trunc(n.yos/2)
    max.NASC.b = max.NASC.m/5 # maximum bin depth ensonified
    NASC.dat <- NASC.leg1[1:max.NASC.b,match(yo.locs,colnames(NASC.leg1))]
    
    # NASC.dat is the field of NASC values the glider is sampling from
    colnames(NASC.dat) <- paste('X',colnames(NASC.dat),sep='')
    ship.Sv[iyr] <- mean(apply(NASC.dat,2,sum,na.rm=TRUE))

    d.bins <- length(unique(out$z.bin)) # depth bins in glider path
    h.bins <- length(unique(out$d.bin)) # horizontal bins
    smpl.pth <- array(dim=list(d.bins,h.bins))
		                # smpl.pth is the glider sampling matrix
                                # calculated in gldry.r.
                                # It is the same for each replicate
                                # but starts in different columns of NASC.dat
    k <- 0
    for(i in 1:d.bins)
      for(j in 1:h.bins){
        k <- k+1
        smpl.pth[i,j] <- out$dive.id[k]
        }
	  # smpl.pth bins are now populated with either
          # a sample yo number or, for unsampled bins, 'NA'

    # Calculate gldr.dat, the vector of (NASC) bin values that were sampled.
      smpl.pths <- smpl.pth[,1:(ncol(smpl.pth)/2)]
                                                 # each i.gldr sample will be half the total
                                                 # possible yos, with a random starting
                                                 # column in the first half of NASC.dat
      gldr.dat <- array(dim=c((shutoff.azfp.m+azfp.range.m)/bin.depth.m,
                  ncol(smpl.pths),i.gldr,n.rep)) 
                              # gldr.dat is a 4 dimensional array,
                              # depth(bins) X horizontal dist X passes X replicates
      c.gldrs <- array(dim=c((shutoff.azfp.m+azfp.range.m)/bin.depth.m,ncol(gldr.dat)*i.gldr,n.rep))
      # c.gldrs is a 3 dimensional array of the depths X areas sampled by all gliders X replicates
      gldr.table <- array(dim=c(n.rep,length(qntl.vals)+2))
      dimnames(gldr.table) <- list(1:n.rep,c(qntl.vals,'sum','sd'))
      y.means <- array(dim=c(nrow(gldr.dat),yo.count[iyr],n.rep)) # mean gldr densities at depth
                                                                  # for each yo
      yo.sums <- array(dim=c(yo.count[iyr],n.rep)) # summed depth means for each yo

      for(i.rep in 1:n.rep){
        #set.seed()
        if(smpl.st > 0) { # begin smpl.st > 0 (assign glider starting locations from file)
            gldr.strt.row.no <-  
               read.table(paste('gldr_strt/',i.gldr,'_', AMLR.area,'_',NASC.yrs[iyr],
               '_',depths[i.depth],'_gldr_strt_nos.txt',sep=''),header = TRUE)
            gldr.strt.row.name <- colnames(NASC.dat)[unlist(gldr.strt.row.no)]
            dim(gldr.strt.row.name) <- dim(gldr.strt.row.no)
            gldr.strt.table[,i.rep] <- unlist(gldr.strt.row.no[i.rep]) 
            gldr.strt.table.nm[,i.rep] <- unlist(gldr.strt.row.name[i.rep]) 
        
            #for(i.gldr in n.gldr){
            #   if glider is not going below 500 m (NASC data limit):
            if((shutoff.azfp.m+azfp.range.m)/bin.depth.m <= nrow(NASC.dat)){ 
              gldr.dat[1:nrow(gldr.dat),1:ncol(smpl.pths),i.gldr,i.rep] <- 
                NASC.dat[1:nrow(gldr.dat),as.numeric(gldr.strt.row.no[i.gldr,i.rep]):
                (as.numeric(gldr.strt.row.no[i.gldr,i.rep])+ncol(smpl.pths)-1)]
              }
            #   if glider is going below 500 m:
            if((shutoff.azfp.m+azfp.range.m)/bin.depth.m > nrow(NASC.dat)){
              gldr.dat[1:nrow(NASC.dat),1:ncol(smpl.pths),i.gldr,i.rep] <- 
                NASC.dat[1:nrow(NASC.dat),as.numeric(gldr.strt.row.no[i.gldr,i.rep]):
                (as.numeric(gldr.strt.row.no[i.gldr,i.rep])+ncol(smpl.pths)-1)]
              }
            #} # end of inside i.gldr loop
           }else {  # end smpl.st > 0,begin smpl.st == 0 (random glider starting positions)
             gldr.strt.row.no <- sample(1:(ncol(NASC.dat)/2),i.gldr,replace=TRUE)
             gldr.strt.row.name <- colnames(NASC.dat)[unlist(gldr.strt.row.no)]
             dim(gldr.strt.row.name) <- dim(gldr.strt.row.no)
             # starting position is randomized for each n.gldr for the replicate in line above

             if(i.gldr>1){ 
             gldr.strt.table[,i.rep] <- gldr.strt.row.no 
             gldr.strt.table.nm[,i.rep] <- gldr.strt.row.name 
             } 
             else { # i.gldr==1
             gldr.strt.table[i.gldr,i.rep] <- gldr.strt.row.no 
             gldr.strt.row.name <- colnames(NASC.dat)[unlist(gldr.strt.row.no)]
             dim(gldr.strt.row.name) <- dim(gldr.strt.row.no)
             gldr.strt.table.nm[,i.rep] <- gldr.strt.row.name 
             }
            #for(i.gldr in 1:i.gldr){
             #   if glider is not going below 500 m (NASC data limit):
            if((shutoff.azfp.m+azfp.range.m)/bin.depth.m <= nrow(NASC.dat)){ 
            for (g in 1:i.gldr)
             gldr.dat[1:nrow(gldr.dat),1:ncol(smpl.pths),g,i.rep] <- 
               NASC.dat[1:nrow(gldr.dat),gldr.strt.row.no[g]:
                       (gldr.strt.row.no[g]+ncol(smpl.pths)-1)]
             }
            #   if glider is going below 500 m:
            if((shutoff.azfp.m+azfp.range.m)/bin.depth.m > nrow(NASC.dat)){
            for (g in 1:i.gldr)
             gldr.dat[1:nrow(NASC.dat),1:ncol(smpl.pths),g,i.rep] <- 
               NASC.dat[1:nrow(NASC.dat),gldr.strt.row.no[g]:
                       (gldr.strt.row.no[g]+ncol(smpl.pths)-1)]
             }
           # } # end of i.gldr loop        
           } # end smpl.st == 0

      # combine samples from multiple gliders per replicate
      tmp.dat <- cbind(gldr.dat[1:nrow(gldr.dat),1:ncol(gldr.dat),,i.rep])
      dim(tmp.dat) <- c(nrow(gldr.dat),ncol(gldr.dat)*i.gldr)
      c.gldrs[,,i.rep] <- tmp.dat # c.gldrs is the combined glider samples for the replicate
      smpl.pths.multi <- smpl.pths # create sampling mask for multiple gliders
      if(i.gldr > 1)
        for(j in 2:i.gldr)
          smpl.pths.multi <- cbind(smpl.pths.multi,smpl.pths)

      for(i.yo in 1:max(yo.count[iyr],na.rm=TRUE))
        for(d.bin in 1:nrow(gldr.dat)){
          y.means[d.bin,i.yo,i.rep] <- mean(c.gldrs[d.bin,,i.rep]
                                [which(smpl.pths.multi[d.bin,]==i.yo)],na.rm=TRUE)
          }
      yo.sums <- apply(y.means,c(2,3),sum,na.rm=TRUE) # summed depth means for each yo
      gldr.sums[i.rep,iyr] <- mean(yo.sums[,i.rep],na.rm=TRUE) # mean of yo sums for the year
      gldr.sds[i.rep,iyr] <- sd(yo.sums[,i.rep],na.rm=TRUE)
      gldr.table[i.rep,] <- c(quantile(c.gldrs[,,i.rep][smpl.pths>0],
        probs = qntl.vals,na.rm=TRUE),mean(yo.sums[,i.rep],na.rm=TRUE),
        sd(yo.sums[,i.rep],na.rm=TRUE)
        )
      # gldr.table contains the quantiles and mean and sd of yo sums
      # for the combined gliders in each replicate      
      #####################
      write.table(smpl.pths.multi,paste('tables/smpl_paths/',i.gldr,'_',
                                 AMLR.area,'_',smpl.yr,'_',depths[i.depth],'_spths.txt',sep='')) 
      } # end of i.rep loop

    write.table(round(c.gldrs,5),paste('tables/cgldrs/',i.gldr,'_',
                                 AMLR.area,'_',smpl.yr,'_',depths[i.depth],'_cgldrs.txt',sep='')) 

     # Calculate quantiles for glider and ship NASC values for the year
    NASC.table <- c(quantile(NASC.dat,
      probs = qntl.vals,na.rm=TRUE),mean(apply(NASC.dat,2,sum,na.rm=TRUE)),sd(apply(NASC.dat,2,sum,na.rm=TRUE)))
    ship.table[iyr,] <- c(NASC.yrs[iyr],NASC.table)

    if(save.tables == 1){
      write.table(round(gldr.table,digits=1),
        file=paste('tables/gldrtables/',i.gldr,'_',
         AMLR.area,'_',smpl.yr,'_',depths[i.depth],'m_gldrtable.txt',sep=''))
      write.table(as.matrix(gldr.strt.table.nm),paste('tables/gldr_strt/',
         i.gldr,'_', AMLR.area,'_',NASC.yrs[iyr],'_',depths[i.depth],'_gldr_strt_names.txt',sep=""))
      write.table(gldr.strt.table,paste('tables/gldr_strt/',
         i.gldr,'_', AMLR.area,'_',NASC.yrs[iyr],'_',depths[i.depth],'_gldr_strt_nos.txt',sep=''))
      write.table(round(yo.sums,5),paste('tables/yo_sums/',
         i.gldr,'_', AMLR.area,'_',smpl.yr,'_',depths[i.depth],'_yo_sums.txt',sep=''))
      write.csv(round(t(NASC.dat),4),
         file=paste('tables/NASC_yrs/',i.gldr,'_',AMLR.area,'_',smpl.yr,'_NASC_dat.csv',sep=''))
      }
      #####################
    } # end iyr loop
  ship.table[,2:8] <- round(ship.table[,2:8],digits=1)

  if(save.tables == 1){
    write(ship.Sv,paste('tables/ship_Sv_',AMLR.area,'.txt'),ncolumns=10)
    write.table(ship.table,row.names=FALSE,
        file=paste('tables/shiptable',AMLR.area,'.txt',sep=''))
    write.table(round(gldr.sums,5),paste('tables/gldr_stats/',
        i.gldr,'_', AMLR.area,'_',depths[i.depth],'_gldr_sums.txt',sep=''))
    write.table(round(gldr.sds,5),paste('tables/gldr_stats/',
        i.gldr,'_', AMLR.area,'_',depths[i.depth],'_gldr_sd.txt',sep=''))
    write(yo.count*i.gldr,paste('tables/yo_count/',
        i.gldr,'_', AMLR.area,'_',depths[i.depth],'_yo_count.txt',sep=''),ncolumns=10)
  } # end save.tables
 } # end of i.depth, n.gldr loops
} # end function gldrs()
