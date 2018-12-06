#svm
library(kernlab)
set.seed(64532)

#caso separable
n1 <- 50;mx1 <- -1;my1 <- 1;s1=0.5
n2 <- 50;mx2 <- 1;my2 <- -1;s2=0.5
a1 <- matrix(c(rnorm(n1,mx1,s1),rnorm(n1,my1,s1)),ncol=2)
a2 <- matrix(c(rnorm(n2,mx2,s2),rnorm(n2,my2,s2)),ncol=2)
dat <- rbind(a1,a2)
rx <- range(dat[,1]); rx <- rx+c(-1,1)*0.05*(rx[2]-rx[1])
ry <- range(dat[,2]); ry <- ry+c(-1,1)*0.05*(ry[2]-ry[1])
plot(0,0,type="n",xlim = rx,ylim = ry,xaxt="n",yaxt="n",xlab="",
     ylab="",main="Hiperplano Separador",cex.main=1)
points(a1,col="blue",pch=19)
points(a2,col="red",pch=19)
yy <- c(rep(1,n1),rep(-1,n2))
d <- 3
n <- n1+n2
M <- 10000
cc <- rep(0,d)
HH <- matrix(0,d,d)
HH[2,2] <- HH[3,3] <- 1
#HH=diag(d)
AA <- cbind(yy,yy*dat)
bb <- rep(1,n)
ll <- rep(-M,d)
uu <- rep(M,d)
rr <- rep(M,n)
out <- ipop(c=cc,H=HH,A=AA,b=bb,l=ll,u=uu,r=rr)
bas <- slot(out,"primal")
xs <- seq(rx[1],rx[2],length=10)
ys <- -(bas[1]+xs*bas[2])/bas[3]
desp <- sqrt(1+(bas[2]/bas[3])^2)/sqrt(bas[2]^2+bas[3]^3)
desp <- 1/abs(bas[3])
ysu <- ys+desp
ysd <- ys-desp
lines(xs,ys,lwd=2)
lines(xs,ysu,lty=2)
lines(xs,ysd,lty=2)
#vectores soporte
lamb <- slot(out,"dual")
sop <- which(lamb>0.00000001)
points(dat[sop,],pch=16)

#############
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
plot(0,0,type="n",xlim=ry,ylim=ry,xaxt="n",yaxt="n",xlab="",ylab="")
points(aa,col="blue",pch=24,bg="blue",cex=0.8)
points(bb,col="blue",pch=24,bg="blue",cex=0.8)
points(cc,col="red",pch=19,cex=0.8)
points(dd,col="red",pch=19,cex=0.8)
points(ee,col="red",pch=19,cex=0.8)
abline(h=0,v=0,col=gray(0.9),cex=0.8)
####SVM primal sin transformaciones
CC <- 4
hh <- CC*c(rep(0,3),rep(1,n))
HH <- diag(c(0,1,1,rep(0,n)))
AA <- cbind(y*cbind(rep(1,n),dat),diag(n))
b <- rep(1,n)
M <- 1000
r <- rep(M,n)
l <- c(rep(-M,3),rep(0,n))
u <- rep(M,3+n)
outk <- ipop(hh,HH,AA,b,l,u,r)
bas <- (primal(outk))[1:3]
xs <- seq(rx[1],rx[2],length=2)
ys <- -(bas[1]+xs*bas[2])/bas[3]
desp <- 1/abs(bas[3])
ysu <- ys+desp
ysd <- ys-desp
lines(xs,ys,lwd=2,col="black")
lines(xs,ysu,lty=2)
lines(xs,ysd,lty=2)
alf <- dual(outk)
ww <- which((alf>0.01)&(alf<3.99))
points(dat[ww,1],dat[ww,2],pch=16)
bet <- colSums(alf*y*dat)
b0 <- mean(y[ww]-dat[ww,]%*%bas[2:3])
#resolviendo el problema dual
CC <- 4
hh <- rep(-1,n)
HH <- (y*dat)%*%t(y*dat)
AA <- t(as.matrix(y))
b <- 0
r <- 0
l <- rep(0,n)
u <- rep(CC,n)
out <- ipop(hh,HH,AA,b,l,u,r)
alf <- primal(out)
bet <- colSums(alf*y*dat)
ww <- which((alf>0.01)&(alf<3.99))
b0 <- mean(y[ww]-dat[ww,]%*%bas[2:3])

#a que distancia estan esos del hiperplano
(b0+dat[ww,]%*%bet)/sqrt(bet[1]^2+bet[2]^2)

#margen
marg <- 1/sqrt(bet[1]^2+bet[2]^2)