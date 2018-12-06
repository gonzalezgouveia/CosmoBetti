#SVM con kernel para datos de esferas y cubos
library(kernlab)
library(TDA)
n <- 50 #cantidad de puntos de la esfera
B <- 500 #cantidad de nubes de datos. Solo poner cantidades pares
dat1 <- numeric(0)
for(i in 1:(B/2)){
  a <- sphereUnif(n,2)
  dat1 <- rbind(dat1,matrix(a,nrow=1))
}
dat2 <- numeric(0)
for(i in 1:(B/2)){
  #a <- torusUnif(n,1,3)
  #a <- rcubo(n,r=(pi/6)^(1/3))
  a <- sphereUnif(n,2)
  dat2 <- rbind(dat2,matrix(a,nrow=1))
}
dat <- rbind(dat1,dat2)
y <- c(rep(-1,B/2),rep(1,B/2))
d <- n*3 #150
n <- B
#empiezan los calculos
#SVM con kernel base radial
bet <- 0.05 #parametros bet y CC de validacion cruzada
CC <- 1
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

#matriz de confusion
#prediccion bajo el modelo
VV2 <- kernelPol(rbf,x=dat,y=dat,z=aas*y,k=rep(1,n))
rr <- sign(colSums(VV2)+tet0)
conf11 <- length(which(rr[1:(n/2)]>0))  #bien clasificados del grupo 1
conf12 <- length(which(rr[1:(n/2)]<0))  #mal clasificados del grupo 1
conf21 <- length(which(rr[(n/2+1):n]>0))  #mal clasificados del grupo 2
conf22 <- length(which(rr[(n/2+1):n]<0))  #bien clasificados del grupo 2
confmat <- matrix(c(conf11,conf12,conf21,conf22),ncol=2,byrow=T)
print(confmat)
errLDA <- 100*(1-sum(diag(confmat))/sum(confmat))
print(errLDA)


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

#validacion cruzada 10-fold cv
aux <- sample(1:n)
kut <- seq(0,n,by=(n/10))
bet <- c(.01,.05)
CC <- c(1,10,100)
nbet <- length(bet);nCC <- length(CC); CVerr <- matrix(0,nbet,nCC)
datos <- cbind(y,dat)
for(ii in 1:nbet){
  print(paste("ii=",ii))
  for(jj in 1:nCC){
    print(paste("jj=",jj))
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
auxind <- which(CVerr == min(CVerr),arr.ind = T)[1,]
print(auxind)
bet[auxind[1]]
CC[auxind[2]]
