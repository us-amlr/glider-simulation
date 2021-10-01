# for Fig. 2a: yo.depths.m=rep(250,3); max.z.m=255
# for Fig. 2b: yo.depths.m=rep(1000,3); max.z.m=1005

library(lattice)
Fig2 <- function(yo.depths.m=rep(250,3), # max yo depth, no.yos
angle.deg=22.7,
azfp.range.m=100,
shutoff.azfp.m=150,
bin.length.m=100,
bin.depth.m=5,
max.z.m=255,  #  this value sets the bottom of the plot (255 or 1005)
ship.effective.m=250){
  # convert degrees to radians for use in trig functions
  angle.rad<-angle.deg*pi/180
  # z.rate is the rate at which the glider dives
  # z.rate is in number of bins through which dive progresses per horizontal bin traversed
  z.rate<-ceiling(tan(angle.rad)*bin.length.m/bin.depth.m)

  n.yos<-length(yo.depths.m)  
  # set up the bins (all distances in numbers of bins)
  yo.dist<-rep(NA,n.yos)
  yo.dist.start<-c(1,rep(NA,n.yos-1))
  yo.dist.mid<-rep(NA,n.yos)
  yo.dist.end<-rep(NA,n.yos)
  
  for(i in 1:n.yos){
    yo.dist[i]<-ceiling(2*(yo.depths.m[i]/bin.depth.m)/z.rate)
    yo.dist.end[i]<-sum(yo.dist[1:i])
    if(i>=2){
      yo.dist.start[i]<-yo.dist.end[i-1]+1
    }
    yo.dist.mid[i]<-ceiling((yo.dist.end[i]+yo.dist.start[i])/2)
  }
  
  d.seq<-1:sum(yo.dist)
  z.seq <- 1:(max.z.m/bin.depth.m)
  
  out<-expand.grid(d.bin=d.seq, z.bin=z.seq)
  out$dive.id<-rep(NA,dim(out)[1])
  
  # identify which bins are ensonified during each dive (descending leg of half-yo)
  for(i in 1:n.yos){
    dist.range<-yo.dist.start[i]:yo.dist.mid[i]
    z.start<-rep(NA,length(dist.range))
    z.end<-rep(NA,length(dist.range))
    z.ext<-rep(NA,length(dist.range))
    
    azfp.limit<-min(yo.depths.m[i]/bin.depth.m,shutoff.azfp.m/bin.depth.m)
    
    break.me<-FALSE
    for(j in 1:length(dist.range)){
      if(j==1){
        z.start[j]<-1
        
        if(z.rate<azfp.limit){
          z.end[j]<-z.rate
        } else {
          z.end[j]<-azfp.limit
          break.me<-TRUE
        }        
      }else{
        z.start[j]<-z.end[j-1]+1
        
        if((z.start[j]+z.rate)<azfp.limit){
          z.end[j]<-z.start[j]+z.rate  
        } else {
          z.end[j]<-azfp.limit
          break.me<-TRUE
        }        
      }
      
      tt<-(out$z.bin%in%(z.start[j]:z.end[j]))&(out$d.bin==(yo.dist.start[i]+j-1))
      out$dive.id[tt]<-i
      
      z.ext[j]<-z.end[j]+(azfp.range.m/bin.depth.m)  
      ttt<-(out$z.bin%in%((z.end[j]+1):z.ext[j]))&(out$d.bin==(yo.dist.start[i]+j-1))
        out$dive.id[ttt]<-n.yos+1  
        
      if(break.me){break}
    } # end j loop
  } # end i loop 

# set up the color scheme
    cols<- c(rep('black',3),'light gray')
    cols.at<-0:(n.yos+1)
  
  Xlab<-'Surface Distance (m)'
  Ylab<-'Depth (m)'
  
out.plot<-levelplot(dive.id~d.bin*100+I(-z.bin)*5,data=out,col.regions=cols,at=cols.at,
 		      scales=list(x=list(cex=0.9),y=list(cex=0.9)),
		      colorkey=FALSE,
                      xlab=list(Xlab,cex=0.9),
                      ylab=list(Ylab,cex=0.9))
 print(out.plot)

} # end function(Fig2)
