gldry <-function(std.yo.depth.m=200,shutoff.azfp.m=200,max.z.m=200,
        yo.locs=rep('a',25)){
angle.deg=22.7
# std.yo.depth.m=200 # assigned in calling script
azfp.range.m=100
#shutoff.azfp.m=200 # assigned in calling script
view.azfp.range=FALSE
bin.length.m=100
bin.depth.m=5
# max.z.m=200 # assigned in calling script
max.z.b <- max.z.m/bin.depth.m
ship.effective.m=250
  # convert degrees to radians for use in trig functions
  angle.rad<-angle.deg*pi/180
  # z.rate is the rate at which the glider dives
  # z.rate is in number of vertical bins through which dive progresses per horizontal bin traversed
  z.rate<-ceiling(tan(angle.rad)*bin.length.m/bin.depth.m) # 9 depth bins per Interval
  #
  #
  #n.yos <- ceiling(length(yo.locs)/ceiling(max.z.b/z.rate))
  n.yos <- (length(yo.locs)*bin.length.m)/(ceiling(2*(std.yo.depth.m/bin.depth.m)/z.rate)*100)
  #yo.depths.m=rep(max.z.m,n.yos)
  yo.depths.m=rep(std.yo.depth.m,n.yos)

  # n.yos<-length(yo.depths.m)  
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
    # when view.azfp.range=FALSE ensonified bins are identified by integer values for dive.id
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
      if(view.azfp.range){
        out$dive.id[ttt]<-n.yos+1  
      } else {
        out$dive.id[ttt]<-i
      }
        
      if(break.me){break}
    } # end j loop
  } # end i loop 
return(list(n.yos,out,azfp.range.m,bin.depth.m))
} # end gldry()