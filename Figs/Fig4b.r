# Fig 4b
DF <- data.frame(NASC.leg1[1:50,])
prop.0 <- apply(DF == 0, 1, sum,na.rm=TRUE)/(
          apply(DF == 0, 1, sum,na.rm=TRUE)+apply(DF > 0, 1, sum,na.rm=TRUE))

depths <- rev(c(1:50))

par(cex=0.9,cex.lab=0.9)
plot(x=rev(prop.0),y=1:length(depths),axes=FALSE,
     pch=19,ylab= "Depth (m)",
     xlab= expression(paste("Proportion of zero ",s[a]," at depth"))
    )
axis(1)
axis(2,at = length(depths):1,labels=rev(depths)*5)
box()
