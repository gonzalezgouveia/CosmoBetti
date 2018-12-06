#caso no separable
library(kernlab)
set.seed(656431)
#par(mfrow=c(2,1),mar=c(2,1,1,1))
n1 <- 50;mx1 <- -1; my1 <-  1;s1=0.6
n2 <- 10;mx2 <-  0; my2 <-  0;s2=0.6
n3 <- 20;mx3 <- -1; my3 <- -1;s3=0.6
n4 <- 20;mx4 <-  1; my4 <- -1;s4=0.6
n5 <- 20;mx5 <-  1; my5 <-  1;s5=0.6
n <- n1+n2+n3+n4+n5
y <- c(rep(1,n1),rep(1,n2),rep(-1,n3),rep(-1,n4),rep(-1,n5))
d <- 2
aa <- matrix(c(rnorm(n1,mx1,s1),rnorm(n1,my1,s1)),ncol=d)
bb <- matrix(c(rnorm(n2,mx2,s2),rnorm(n2,my2,s2)),ncol=d)
cc <- matrix(c(rnorm(n3,mx3,s3),rnorm(n3,my3,s3)),ncol=d)
dd <- matrix(c(rnorm(n4,mx4,s4),rnorm(n4,my4,s4)),ncol=d)
ee <- matrix(c(rnorm(n5,mx5,s5),rnorm(n5,my5,s5)),ncol=d)
dat <- rbind(aa,bb,cc,dd,ee)
rx <- range(c(aa[,1],bb[,1],cc[,1],dd[,1],ee[,1]))
ry <- range(c(aa[,2],bb[,2],cc[,2],dd[,2],ee[,2]))
#SVM con kernel base radial
bet <- 0.4 #parametros bet y CC de validacion cruzada
CC <- 3
rbf <- rbfdot(sigma=bet/2)
HH <- kernelPol(rbf,dat,z=y)
gg <- -rep(1,n)
AA <- t(y)
b <- 0
r <- 0
l <- rep(0,n)
u <- rep(CC,n)
outk <- ipop(gg,HH,AA,b,l,u,r)
aas <- primal(outk)
indx <- which(aas>0.001 & aas<0.999*CC)
nm <- length(indx)
UU <- kernelPol(rbf,x=dat,y=dat[indx,],z=aas*y,k=rep(1,nm))
tet0 <- median(y[indx]-colSums(UU))
M <- 100
xx <- seq(rx[1],rx[2],length=M)
yy <- seq(ry[1],ry[2],length=M)
datt <- as.matrix(expand.grid(xx,yy))
np <- dim(datt)[1]
VV <- kernelPol(rbf,x=dat,y=datt,z=aas*y,k=rep(1,np))
zz <- colSums(VV)+tet0
zz <- matrix(zz,ncol=M)
plot(0,0,type="n",xlim=ry,ylim=ry,xaxt="n",yaxt="n",xlab="",ylab="")
image(xx,yy,zz,levels=0,add=T,drawlabels=F,lwd=2)
contour(xx,yy,zz,levels=0,add=T,drawlabels = F,lwd=2)
contour(xx,yy,zz,levels=c(-1,1),add=T,drawlabels = F,lwd=1)
points(aa,col="blue",pch=24,bg="blue",cex=0.8)
points(bb,col="blue",pch=24,bg="blue",cex=0.8)
points(cc,col="red",pch=19,cex=0.8)
points(dd,col="red",pch=19,cex=0.8)
points(ee,col="red",pch=19,cex=0.8)
points(dat[indx,],col="black",pch="o") #vectores soporte

#los margenes estan en -1 y 1, como se ve de los siguientes calculos:
np <- length(indx)
datt <- dat[indx,]
VV <- kernelPol(rbf,x=dat,y=datt,z=aas*y,k=rep(1,np))
rr <- colSums(VV)+tet0

#matriz de confusion
#prediccion bajo el modelo
VV2 <- kernelPol(rbf,x=dat,y=dat,z=aas*y,k=rep(1,n))
rr <- sign(colSums(VV2)+tet0)
conf11 <- length(which(rr[1:(n/2)]>0))  #bien clasificados del grupo 1
conf12 <- length(which(rr[1:(n/2)]<0))  #mal clasificados del grupo 1
conf21 <- length(which(rr[(n/2):n]>0))  #mal clasificados del grupo 2
conf22 <- length(which(rr[(n/2):n]<0))  #bien clasificados del grupo 2
confmat <- matrix(c(conf11,conf12,conf21,conf22),ncol=2,byrow=T)
print(confmat)
errLDA <- 100*(1-sum(diag(confmat))/sum(confmat))
print(errLDA)

#validacion cruzada. seleccion de beta y C, linea 22 y 23
#funcion CV
CV <- function(datos,sel,bet,CC){
  dat <- datos[sel,-1]
  y <- datos[sel,1]
  n <- dim(dat)[1]
  rbf <- rbfdot(sigma =bet/2)
  HH <- kernelPol(rbf,dat,z=y)
  gg <- -rep(1,n)
  AA <- t(y)
  b <- 0
  r <- 0
  l <- rep(0,n)
  u <- rep(CC,n)
  outk <- ipop(gg,HH,AA,b,l,u,r)
  aas <- primal(outk)
  indx <- which(aas>0.001 & aas<0.999*CC)
  nm <- length(indx)
  UU <- kernelPol(rbf,x=dat,y=dat[indx,],z=aas*y,k=rep(1,nm))
  tet0 <- median(y[indx]-colSums(UU)) #calculo de teta0
  dt <- datos[-sel,-1]
  ndt <- dim(dt)[1]
  VV <- kernelPol(rbf,x=dat,y=dt,z=aas*y,k=rep(1,ndt))
  zz <- colSums(VV)+tet0
  prd <- sign(zz)
  return(prd)
}

set.seed(65656)
aux <- sample(1:n)
kut <- seq(0,n,by=12)
bet <- c(.05,.1,.2,.3,.4,.5,.7,1,3,8,20)
CC <- c(1,2,3,4,5,10,100)
nbet <- length(bet);nCC <- length(CC); CVerr <- matrix(0,nbet,nCC)
datos <- cbind(y,dat)
for(ii in 1:nbet){
  for(jj in 1:nCC){
    err <- rep(0,10)
    for(i in 1:10){
      sel <- aux[-((kut[i]+1):(kut[i+1]))]
      pre <- CV(datos,sel,bet[ii],CC[jj])
      err[i] <- sum((y[-sel]*pre)==-1)/12
    }
    CVerr[ii,jj] <- 100*mean(err)
  }
}
print(round(CVerr,3))
