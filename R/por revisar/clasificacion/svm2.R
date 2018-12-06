#simulando los datos
n <- 50 #cantidad de puntos de la esfera
B <- 500 #cantidad de nubes de datos. Solo poner cantidades pares
dat1 <- numeric(0)
for(i in 1:(B/2)){
  a <- sphereUnif(n,2)
  #a <- rcubo(n,r=(pi/6)^(1/3))
  dat1 <- rbind(dat1,matrix(a,nrow=1))
}
dat2 <- numeric(0)
for(i in 1:(B/2)){
  #a <- torusUnif(n,1,3)
  #a <- rcubo(n,r=(pi/6)^(1/3))
  #a <- sphereUnif(n,2)
  a <- klein.unif(k=n,sigma=0.0)
  for(j in 1:3){
    a[,j] <- a[,j]/(max(a[,j])-min(a[,j]))*2
  }
  dat2 <- rbind(dat2,matrix(a,nrow=1))
}
dat <- rbind(dat1,dat2)
y <- c(rep(0,B/2),rep(1,B/2))
#mi version del svm que imprime las matrices de confusion
library(kernlab)
#n <- 500
#y <- c(rep(1,n/2),rep(-1,n/2))
#d <- 150
#aa <- matrix(rnorm(d*n/2,1,1),ncol=d)
#bb <- matrix(rnorm(d*n/2,0.8,1),ncol=d)
#dat <- rbind(aa,bb)
#despues de correr el de clasi2.R que tiene lo de LDA
d <- n*3 #150
n <- B
y[1:(B/2)] <- -1
####SVM primal sin transformaciones
CC <- 4
hh <- CC*c(rep(0,(d+1)),rep(1,n))
HH <- diag(c(0,rep(1,d),rep(0,n)))
AA <- cbind(y*cbind(rep(1,n),dat),diag(n))
b <- rep(1,n)
M <- 1000
r <- rep(M,n)
l <- c(rep(-M,(d+1)),rep(0,n))
u <- rep(M,(d+1)+n)
outk <- ipop(hh,HH,AA,b,l,u,r)#solves quadratic programming problem
print(how(outk))
bas <- primal(outk)[1:((d+1))]
#buscando los soportes
alf <- dual(outk)
ww <- which((alf>0.01)&(alf<3.99))
#points(dat[ww,1],dat[ww,2],pch=16)
bet <- colSums(alf*y*dat)
b0 <- mean(y[ww]-dat[ww,]%*%bas[2:(d+1)])

#como saber si un punto en dat esta clasificado por uno o por otro
#los que son igual a cero es porque estan mal clasificados
#su clasificacion era 1 y la prediccion dio -1 o viceversa
#los que dieron 2 o -2 es porque estan bien clasificados
#su clasificacion era 1 y la preduccuon dio 1 o viceversa
conf11 <- sum(y[1:(n/2)]+sign(b0+dat[1:(n/2),]%*%bet)!=0) #bien clasificados del grupo 1
conf12 <- sum(y[1:(n/2)]+sign(b0+dat[1:(n/2),]%*%bet)==0) #mal clasificados del grupo 1
conf21 <- sum(y[(n/2+1):n]+sign(b0+dat[(n/2+1):n,]%*%bet)==0) #bueno clasificados del grupo 1
conf22 <- sum(y[(n/2+1):n]+sign(b0+dat[(n/2+1):n,]%*%bet)!=0) #bueno clasificados del grupo 1
confmat <- matrix(c(conf11,conf12,conf21,conf22),ncol=2,byrow=T)
print(confmat)
errLDA <- 100*(1-sum(diag(confmat))/sum(confmat))
print(errLDA)