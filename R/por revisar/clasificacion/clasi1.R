library(ElemStatLearn)
library(mvtnorm)
data("vowel.test")
data("vowel.train")
head(vowel.train,15)
dat <- vowel.train
y <- dat[,1] #las clasificaciones
dat <- dat[,-1]
d <- 10#cantidad de variables
c <- 11#cantidad de clases
#LDA
medias <- by(dat,y,colMeans)
vars <- by(dat,y,var)
ns <- by(dat,y,function(x) dim(x)[1])
Sp <- matrix(0,d,d)
for(i in 1:c){Sp <- Sp+(ns[[i]]-1)*vars[[i]]} #EL ONCE 11
Sp <- Sp/(sum(ns)-c)  #EL ONCE 11
qs <- matrix(0,dim(dat)[1],c)
for(j in 1:c){qs[,j] <- dmvnorm(dat,mean=medias[[j]],sigma=Sp)}
mm <- apply(qs,1,max)
clasepred <- rep(0,dim(dat)[1])
for(i in 1:dim(dat)[1]){
  aa <- which(qs[i,]==mm[i])
  nn <- length(aa)
  if(nn==1){clasepred[i]=aa
  }else{clasepred[i]=aa[sample(1:nn,size=1)]}
}
aa = table(y,clasepred)
errLDA <- 100*(1-sum(diag(aa))/sum(aa))

#EVUALUACION FDA
datP <- vowel.test
yP <- datP[,1]
datP <- datP[,-1]
qsP <- matrix(0,462,11)
for(j in 1:11){qsP[,j] <- dmvnorm(datP,mean=medias[[j]],sigma=Sp)}
mmP <- apply(qsP,1,max)
clasepreP <- rep(0,462)
for(i in 1:462){
  aa <- which(qsP[i,]==mmP[i]) 
  nn <- length(aa)
  if(nn==1){clasepreP[i]=aa
  }else{clasepreP[i]=aa[sample(1:nn,size=1)]}
}
aa <- table(yP,clasepreP)  
aa
errLDAP <- 100*(1-sum(diag(aa))/sum(aa))

