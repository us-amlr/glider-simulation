par(cex=0.9,cex.lab=0.9)
plot(rev(apply(data.frame(NASC.leg1[1:50,]),1,sd,na.rm=T)/
       apply(data.frame(NASC.leg1[1:50,]),1,mean,na.rm=T)),1:length(depths),
       axes=FALSE,xlab= expression(paste("CV ",s[a]," at depth")),pch=19,
       ylab= "Depth (m)")
axis(1)
axis(2,at = length(depths):1,labels=rev(depths)*5)
box()
