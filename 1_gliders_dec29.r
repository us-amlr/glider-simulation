# from '2021/1jan/jan5/R_1_test/R_1_gliders_dec29.txt'
# Dec 29 modified sd column in ship.table
# Calculates glider variance using yos a sampling units
# Calls "R_glider_yos.txt" to calculate glider path
#
# changed rndm.reps to n.rep and n.gldr[i.pth] 
# changed "yo.depth" variable to "gldr.dat"
# "c.gldrs" combines gliders
#  c.gldrs[[i.gldr]][[iyr]] <- list()
#  c.gldrs[[i.gldr]][[iyr]][,,i.rep][smpl.pths>0] contains a replicate sample
#  calculates transect (yo) means as t.m,
#  then mean(t.m.means,na.rm=TRUE) for each year
#
# replaced "n.gldr[i.pth][i.gldr]" with "n.gldr[i.pth]"
# replaced "gldr.pths[i.pth][.gldr]" with "gldr.pths[i.pth]"
# added t.sd,
# added year index to t.m and t.sd
# replaced leg 1 with leg1

rm(list=ls())
NASC.yrs <- c(2001:2009,2011) # using newer data
AMLR.area <- "SA"

n.rep <- 500   # number of replicates
n.gldr <- c(1) #,2,3,4,5)  # number of passes through the NASC.dat database per replicate

max.NASC.m = 250  # maximum NASC depth (m) the glider values will be compared to
depths <- c(150,185,200,300,400,500,700,1000) # max yo depths
azfp.off <- c(150,150,150,150,150,150,150,150) # acoustics off depth
max.z.m=250 # previously I was assigning this as maximum yo depth 

# provide the glider depth and azfp shutoff values, and output path
path.in <- "c:/zot/glider/"
path.out <- paste(path.in,"2020/12dec/dec30/",AMLR.area,"/",sep="")
gldr.pths <- c("1gldr/") #,"2gldrs/","3gldrs/","4gldrs/","5gldrs/")
gldr.depths <- c(
                "150m/","185m/","200m/","300m/",
                "400m/","500m/","700m/","1000m/"
                )

qntl.vals <- c(0.97,0.98,0.99,0.999,1)

max.NASC.b <- max.NASC.m/5 # convert NASC depth from meters to bins
# load the database
shipNASC <- paste(path.in,"2019\\8aug\\aug13\\gldr_2001_11.RData",sep="")

load(shipNASC)
ship.Sv <- vector() # ship acoustic mean integrated returns each year

# remove two intervals from Process_ID 753 (see jul17/notes)
#leg1 <- leg1[!(leg1$Process_ID == 753 & leg1$Interval == 117),]
#leg1 <- leg1[!(leg1$Process_ID == 753 & leg1$Interval == 118),]
# remove 9 records with high Sv_max values
leg1 <- leg1[which(leg1[,"Sv_max"] < -30,arr.ind=TRUE),]

#transect.m <- transect.sd <- vector()
#for(i.gldr in n.gldr[i.pth]){
  #c.gldrs[[i.gldr]] <- list()  
for(i.pth in 1:length(gldr.pths))
for(i.depth in 1:length(gldr.depths)){
  path.depth <- paste(path.out,gldr.pths[i.pth],gldr.depths[i.depth],sep="")
  std.yo.depth.m <- depths[i.depth] # max yo depth
  shutoff.azfp.m <- azfp.off[i.depth]

  # assign the area, year, and depths to be sampled
  yo.count <- vector()
  ship.table <- array(dim=c(length(NASC.yrs),length(qntl.vals)+3))
  colnames(ship.table) <- c("year",c(qntl.vals,"sum","sd"))
  options(stringsAsFactors=FALSE)

  gldr.sums <- gldr.sds <- array(dim=c(n.rep,length(NASC.yrs)))
  dimnames(gldr.sums) <- dimnames(gldr.sds) <- list(paste("r.",c(1:n.rep),sep=""),NASC.yrs)
  ship.sums <- vector(length=length(NASC.yrs))
  names(ship.sums)  <- NASC.yrs
  for(iyr in 1:length(NASC.yrs)){
    gldr.strt.table <- array(dim=c(n.gldr[i.pth],n.rep))
    smpl.yr <- c(NASC.yrs[iyr])
    NASC.vals <- leg1[substr(leg1$Date_M,1,4) == smpl.yr &
                    leg1$Depth_mean <= shutoff.azfp.m+100 &
		    leg1$area == AMLR.area,]
    # NASC.vals is the subset of the entire ship database
    # that will be sampled with the glider. It is used to
    # find the Process_ID and Interval values for sampling

    yo.locs <- PID <- vector()  
    # find the process and interval IDs (column names) in NASC.vals
    PID.uniq <- unique(NASC.vals$Process_ID)
    k <- 0
    for (i in 1:length(PID.uniq)){
      PID[i] <- PID.uniq[i]
      PID.int <- unique(NASC.vals$Interval[NASC.vals$Process_ID == PID[i]])
      for(j in 1:length(PID.int)){
        k <- k+1
        yo.locs[k] <- paste(PID.uniq[i],".",PID.int[j],sep="")
        }
      # yo.locs are the column names to be sampled in NASC.leg1
      # NASC.leg1 is the ship NASC matrix from shipNASC.RData
      }
    # use "R_glider_yos.txt" outputs (out$z.bin, out$d.bin, out$dive.id)
    # to calculate the glider locations (smpl.pth) in the NASC grid
    source(paste(path.in,"2020/R_glider_yos.txt",sep=""))
    yo.count[iyr] <- trunc(n.yos/2)
    NASC.dat <- NASC.leg1[1:max.NASC.b,match(yo.locs,colnames(NASC.leg1))]
    # NASC.dat is the field of NASC values the glider is sampling from
    ship.Sv[iyr] <- mean(apply(NASC.dat,2,sum,na.rm=TRUE))

    d.bins <- length(unique(out$z.bin)) # depth bins in glider path
    h.bins <- length(unique(out$d.bin)) # horizontal bins
    smpl.pth <- array(dim=list(d.bins,
                       h.bins))
		                # smpl.pth is the glider sampling matrix
                                # calculated in R_glider_yos.txt.
                                # It is the same for each replicate
                                # but starts in different columns of NASC.dat
    k <- 0
    for(i in 1:d.bins)
      for(j in 1:h.bins){
        k <- k+1
        smpl.pth[i,j] <- out$dive.id[k]
        }
	  # smpl.pth bins are now populated with either
          # a sample yo number or, for unsampled bins, "NA"

    # Calculate gldr.dat, the vector of (NASC) bin values that were sampled.
      smpl.pths <- smpl.pth[,1:(ncol(smpl.pth)/2)]
                                                 # each n.gldr[i.pth] sample will be half the total
                                                 # possible yos, with a random starting
                                                 # column in the first half of NASC.dat
      gldr.dat <- array(dim=c((shutoff.azfp.m+azfp.range.m)/bin.depth.m,
                  ncol(smpl.pths),n.gldr[i.pth],n.rep)) 
                              # gldr.dat is a 4 dimensional array,
                              # depth(bins) X horizontal dist X passes X replicates
      c.gldrs <- array(dim=c((shutoff.azfp.m+azfp.range.m)/bin.depth.m,ncol(gldr.dat)*n.gldr[i.pth],n.rep))
      gldr.table <- array(dim=c(n.rep,length(qntl.vals)+2))
      dimnames(gldr.table) <- list(1:n.rep,c(qntl.vals,"sum","sd"))
      y.means <- array(dim=c(nrow(gldr.dat),yo.count[iyr],n.rep)) # mean gldr densities at depth
                                                                  # for each yo
      yo.sums <- array(dim=c(yo.count[iyr],n.rep)) # summed depth means for each yo

      for(i.rep in 1:n.rep){
        #set.seed(1)
        gldr.strt <- sample(1:(ncol(NASC.dat)/2),n.gldr[i.pth],replace=TRUE) 
        # starting position is randomized for each n.gldr[i.pth] for the replicate in line above
        gldr.strt.table[,i.rep] <- gldr.strt
        for(i.gldr in 1:n.gldr[i.pth]){
     #   if glider is not going below 500 m (NASC data limit):
         if((shutoff.azfp.m+azfp.range.m)/bin.depth.m <= nrow(NASC.dat)){ 
          gldr.dat[1:nrow(gldr.dat),1:ncol(smpl.pths),i.gldr,i.rep] <- 
            NASC.dat[1:nrow(gldr.dat),gldr.strt[i.gldr]:(gldr.strt[i.gldr]+ncol(smpl.pths)-1)]
          }
     #   if glider is going below 500 m:
         if((shutoff.azfp.m+azfp.range.m)/bin.depth.m > nrow(NASC.dat)){
          gldr.dat[1:nrow(NASC.dat),1:ncol(smpl.pths),i.gldr,i.rep] <- 
            NASC.dat[1:nrow(NASC.dat),gldr.strt[i.gldr]:(gldr.strt[i.gldr]+ncol(smpl.pths)-1)]
          }
          } # end of i.gldr loop

      # combine samples from multiple gliders per replicate
      tmp.dat <- cbind(gldr.dat[1:nrow(gldr.dat),1:ncol(gldr.dat),,i.rep])
      dim(tmp.dat) <- c(nrow(gldr.dat),ncol(gldr.dat)*n.gldr[i.pth])
      c.gldrs[,,i.rep] <- tmp.dat # c.gldrs is the combined glider samples for the replicate
      #gldr.table[i.rep,] <- c(quantile(c.gldrs[[i.gldr]][[iyr]][,,i.rep][smpl.pths>0],
      #  probs = qntl.vals,na.rm=TRUE),mean(c.gldrs[[i.gldr]][[iyr]][,,i.rep][smpl.pths>0],na.rm=TRUE),
      #  sd(c.gldrs[[i.gldr]][[iyr]][,,i.rep][smpl.pths>0],na.rm=TRUE))
       # gldr.table contains the quantiles of the combined gliders for each replicate
      for(i.yo in 1:max(yo.count[iyr],na.rm=TRUE))
        for(d.bin in 1:nrow(gldr.dat)){
          y.means[d.bin,i.yo,n.rep] <- mean(c.gldrs[d.bin,,i.rep]
                                [which(smpl.pths[d.bin,]==i.yo)],na.rm=TRUE)
          }
      yo.sums[,i.rep] <- apply(y.means,2,sum,na.rm=TRUE) # summed depth means for each yo
      gldr.sums[i.rep,iyr] <- mean(yo.sums[,i.rep],na.rm=TRUE) # mean of yo sums for the year
      gldr.sds[i.rep,iyr] <- sd(yo.sums[,i.rep],na.rm=TRUE)
      gldr.table[i.rep,] <- c(quantile(c.gldrs[,,i.rep][smpl.pths>0],
        probs = qntl.vals,na.rm=TRUE),mean(yo.sums[,i.rep],na.rm=TRUE),
        sd(yo.sums[,i.rep],na.rm=TRUE)
        )
      # gldr.table contains the quantiles and mean and sd of yo sums
      # for the combined gliders in each replicate      
      #####################
    } # end of i.rep loop

    # Calculate quantiles for glider and ship NASC values for the year
    NASC.table <- c(quantile(NASC.dat,
      probs = qntl.vals,na.rm=TRUE),mean(apply(NASC.dat,2,sum,na.rm=TRUE)),sd(apply(NASC.dat,2,sum,na.rm=TRUE)))
    ship.table[iyr,] <- c(NASC.yrs[iyr],NASC.table)

    write.table(round(gldr.table,digits=1),
        file=paste(path.depth,smpl.yr,"_gldrtable.txt",sep=""))
    write.table(gldr.strt.table,paste(path.depth,"gldr_strt_",NASC.yrs[iyr],".txt",sep=""))
    write.table(round(yo.sums,5),paste(path.depth,smpl.yr,"_",depths[i.depth],"_yo_sums.txt",sep=""))
      #####################
    } # end iyr loop
  write(ship.Sv,paste(path.depth,"_ship_Sv.txt",sep=""),ncolumns=10)

  ship.table[,2:8] <- round(ship.table[,2:8],digits=1)
  write.table(ship.table,row.names=FALSE,
        file=paste(path.depth,"shiptable.txt",sep=""))
  write.table(round(gldr.sums,5),paste(path.depth,depths[i.depth],"_gldr_sums.txt",sep=""))
  write.table(round(gldr.sds,5),paste(path.depth,depths[i.depth],"_gldr_sd.txt",sep=""))
  write(yo.count*n.gldr[i.pth],paste(path.depth,"yo_count.txt",sep=""),ncolumns=10)
  write(gldr.strt,paste(path.depth,"gldr_strt.txt",sep=""),ncolumns=30)
  #write(transect.m,paste(path.depth,"transects_yr_m.txt",sep=""),ncolumns=30)
  #write(transect.sd,paste(path.depth,"transects_yr_sd.txt",sep=""),ncolumns=30)

  } # end of i.depth loop
  #} # end of i.gldr loop
  print(i.depth)

