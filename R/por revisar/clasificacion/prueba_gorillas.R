library(shapes)
library(TDA)
data(apes)
#par(mfrow=c(1,2))

####hace TDA con curvas de Betti a estos datos

plotshapes(apes$x[,,apes$group=="gorm"])
for(i in 1:30){
  points(apes$x[,,apes$group=="gorf"][,,i],col=2)
}
title("Machos negro. Hembras Rojo")
apes$group #hay varias especies de machos y hembras
gorm <- apes$x[,,apes$group=="gorm"]
gorf <- apes$x[,,apes$group=="gorf"]
betti.fun <- function(bcmat,t){
  #return the number of birth - death
  #bcmat es una matriz
  nb <- sum(bcmat[,1]<t) #number of birth
  nd <- sum(bcmat[,2]<t) #number of death
  return(nb-nd)
}
#######curvas de betti
num_figures <- dim(gorm)[3]+dim(gorf)[3]
l <- 500
plm0 <- matrix(0,ncol = l,nrow = num_figures)
plm1 <- matrix(0,ncol = l,nrow = num_figures)
for (pose in 1:num_figures){ # for each pose
  # read the data
  #file <- paste(figures[pose,1],figures[pose,2],sep="")
  print(pose)
  if(pose<=dim(gorm)[3]){
    x <- gorm[,,pose]
  }else{
    x <- gorf[,,(pose-dim(gorm)[3])]
  }
  diag <- alphaShapeDiag(x)$diag
  ph0 <- diag[which(diag[,1]==0),2:3][-1,]
  ph1 <- diag[which(diag[,1]==1),2:3]
  #ph2 <- diag[which(diag[,1]==2),2:3]
  #if(length(ph2)==2){ #cuando hay una sola componente en H2
  #  ph2 <- t(as.matrix(ph2))
  #}
  
  #aqui quiero hacer la cruva de betti para ph0 y para ph1
  xseq <- seq(0,80,length.out = l)
  betticurve0 <- numeric(l)
  betticurve1 <- numeric(l)
  #betticurve2 <- numeric(l)
  for(j in 1:l){
    betticurve0[j] <- betti.fun(ph0,xseq[j])
    betticurve1[j] <- betti.fun(ph1,xseq[j])
    #betticurve2[j] <- betti.fun(ph2,xseq[j])
  }
  plm0[pose,] <- betticurve0
  plm1[pose,] <- betticurve1
  #plm2[pose,] <- betticurve2
}
figures_to_use <- c("macho","hembra")
B = 29
plot(xseq,colMeans(plm0[1:B,]),type="l",main=paste("h0",figures_to_use[1]))
lines(xseq,colMeans(plm0[(B+1):(2*B),]),type="l",col=2)#,main=paste("h0",figures_to_use[2]))
plot(xseq,colMeans(plm1[1:B,]),type="l",main=paste("h1",figures_to_use[1]))
lines(xseq,colMeans(plm1[(B+1):(2*B),]),type="l",col=2)#main=paste("h1",figures_to_use[2]))

