library(ElemStatLearn)
library(mvtnorm)
library(TDA)
n <- 50 #cantidad de puntos de la esfera
B <- 500 #cantidad de nubes de datos. Solo poner cantidades pares
dat1 <- numeric(0)
for(i in 1:(B/2)){
  a <- sphereUnif(n,2)
  #a <- a[order(a[,1]),]
  dat1 <- rbind(dat1,matrix(a,nrow=1))
}
dat2 <- numeric(0)
for(i in 1:(B/2)){
  #a <- torusUnif(n,1,3)
  #a <- rcubo(n,r=(pi/6)^(1/3))
  #a <- sphereUnif(n,2)
  #a <- a[order(a[,1]),]
  a <- klein.unif(k=n,sigma=0.0)
  for(j in 1:3){
    a[,j] <- a[,j]/(max(a[,j])-min(a[,j]))*2
  }
  dat2 <- rbind(dat2,matrix(a,nrow=1))
}
dat <- rbind(dat1,dat2)
y <- c(rep(0,B/2),rep(1,B/2))

#LDA
medias <- by(dat,y,colMeans)
vars <- by(dat,y,var)
ns <- by(dat,y,function(x) dim(x)[1])
Sp <- matrix(0,150,150)
for(i in 1:2){Sp <- Sp+(ns[[i]]-1)*vars[[i]]} #EL ONCE 11
Sp <- Sp/(sum(ns)-2)  #EL ONCE 11
qs <- matrix(0,B,2) #B observaciones 2 clases
for(j in 1:2){qs[,j] <- dmvnorm(dat,mean=medias[[j]],sigma=Sp)}
#sigma=vars[[j]]   #PARA QDA
#sigma=Sp          #PARA LDA
mm <- apply(qs,1,max)
clasepred <- rep(0,B)
for(i in 1:B){
  aa <- which(qs[i,]==mm[i])
  nn <- length(aa)
  if(nn==1){clasepred[i]=aa
  }else{clasepred[i]=aa[sample(1:nn,size=1)]}
}
(aa = table(y,clasepred))
print(aa)
(errLDA <- 100*(1-sum(diag(aa))/sum(aa)))
print(errLDA)
