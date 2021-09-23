# Fig 4a
par(cex=1.4,cex.lab=1.4)
plot(rev(apply(data.frame(NASC.leg1[1:50,]),1,sd,na.rm=T)/
       apply(data.frame(NASC.leg1[1:50,]),1,mean,na.rm=T)),1:length(depths),
       axes=FALSE,xlab= expression(paste("CV ",s[a]," at depth")),pch=19,
       ylab= "Depth (m)")
axis(1)
axis(2,at = length(depths):1,labels=rev(depths)*5)
box()
#dev.new()

# Fig 4b
DF <- data.frame(NASC.leg1[1:50,])
prop.0 <- apply(DF == 0, 1, sum,na.rm=TRUE)/(
          apply(DF == 0, 1, sum,na.rm=TRUE)+apply(DF > 0, 1, sum,na.rm=TRUE))

depths <- rev(c(1:50))

par(cex=1.4,cex.lab=1.4)
plot(x=rev(prop.0),y=1:length(depths),axes=FALSE,
     pch=19,ylab= "Depth (m)",
     xlab= expression(paste("Proportion of zero ",s[a]," at depth"))
    )
axis(1)
axis(2,at = length(depths):1,labels=rev(depths)*5)
box()

