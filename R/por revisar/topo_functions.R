library(spatstat)
library(TDA)

betti.fun <- function(bcmat,t){
  #return the number of birth - death
  #bcmat es una matriz
  nb <- sum(bcmat[,1]<t) #number of birth
  nd <- sum(bcmat[,2]<t) #number of death
  return(nb-nd)
}

meanbc.fun <- function(xmax=0.25,n,l=500,B=1000){
  #meanbc.fun compute the mean betti curve (MBC)
  #for connected components and loops
  #the output is a l=500 times 3 matrix
  #first column is the partition from 0 to xmax
  #second column is the betti curve for connected components
  #third column is the betti curve for loops
  ############3
  #xmax: upper limit of the curve
  #n: number of points of the point process
  #l: size of the output rows for the MBC
  #B: number of repetitios to estimate the MCB
  xseq <- seq(0.0001,xmax^2,length.out = l)
  meanbc0 <- numeric(l)
  meanbc1 <- numeric(l)
  for(i in 1:B){    #repetitions to estimate MCB
    X <- matrix(runif(2*n)*squ,ncol=2)   #binomial process
    diag1 <- alphaComplexDiag(X)$diagram #compute the topological information
    p0  <-  diag1[which(diag1[,1]==0),2:3] #take components
    p1 <- diag1[which(diag1[,1]==1),2:3] #take loops
    for(j in 1:l){   #Update the sum MBC
      meanbc0[j] <- meanbc0[j]+betti.fun(p0,xseq[j])
      meanbc1[j] <- meanbc1[j]+betti.fun(p1,xseq[j])
    }
    #print progress
    if(i%%50==0){print(paste(round(i/B*100,digits = 0),"%"))}
  }
  meanbc0 <- meanbc0/B  #take average
  meanbc1 <- meanbc1/B  #take average
  return(cbind(xseq,meanbc0,meanbc1))
}

distd.fun <- function(meanbc,B = 1000){
  #distd.fun estimate the distribution of distances to the MBC
  l = nrow(meanbc)
  x = meanbc[,1]
  y0 = meanbc[,2]
  y1 = meanbc[,3]
  d0 = numeric(0)      #distance distribution for componentes
  d1 = numeric(0)      #distance distribution for loops
  for(i in 1:B){
    X <- matrix(runif(2*n)*squ,ncol=2)
    diag1 <- alphaComplexDiag(X)$diagram
    p0  <-  diag1[which(diag1[,1]==0),2:3] #components
    p1 <- diag1[which(diag1[,1]==1),2:3] #loops
    z0 = numeric(0)
    z1 = numeric(0)
    daux0 = 0
    daux1 = 0
    for(j in 1:l){
      z0 = betti.fun(p0,x[j]) #number of components at x[j]
      z1 = betti.fun(p1,x[j]) #number of loops at x[j]
      daux0 = daux0 + abs(y0[j]-z0) #distance to MBC components
      daux1 = daux1 + abs(y1[j]-z1) #distance to MBC loops
    }
    d0 = c(d0,daux0) #vector with distribution to distance components
    d1 = c(d1,daux1) #vector with distribution to distance loops
    #print progress
    if(i%%50==0){print(paste(round(i/B*100,digits = 0),"%"))}
  }
  return(cbind(d0,d1))
}
#grafica = plot in spanish
graficabetti <- function(mydata,meanbc,distd,betti,xlim){
  #mydata: a matrix
  #meanbc:  output of meanbc function
  #betti: degree 0 or 1
  #xlim: upper limit of x axis in the plot
  x = meanbc[,1]               #x axis
  y = meanbc[,(betti+2)]           #meanbc_k
  z = numeric(length(x))      #betti curve (not MBC)
  d = 0   #suma de distancias a la media
  n = length(mydata)
  DMD = alphaComplexDiag(mydata)$diagram #topological information
  DMD1 = DMD[which(DMD[,1]==betti),2:3]  #take components or loops
  for(i in 1:length(x)){ #compute the statistic
    z[i] = betti.fun(DMD1,x[i])
    d = d + abs(y[i]-z[i])
  }
  if(d>quantile(distd[,(betti+1)],0.95)){ #testing
    mainhist = "REJECT"
  }else{
    mainhist = "DO NOT REJECT"
  }
  #names for plots
  if(betti == 1){
    main = "Loops. Betti 1"
    ylab = "mean number of loops"
    ylim = 1.5*max(y)
  }else if(betti == 0){
    main = "Components. Betti 0"
    ylab = "mean number of connected components"
    ylim = max(y)
  }
  #plotting
  x11(width =  10,height =  5)
  par(mfrow = c(1,3))
  #x11()
  plot(mydata,main = paste("Binomial Process. n = ",n),
       xlab= paste("number of points = ",n),ylab="")
  
  #x11()
  plot(sqrt(x),y,type = "l",ylab=ylab,col="green",
       main=main,ylim = c(0,ylim),xlim = c(0,xlim))
  for(i in 1:length(x)){
    if(y[i]>z[i]){
      segments(sqrt(x[i]),z[i],sqrt(x[i]),y[i],col="blue")
    }else{
      segments(sqrt(x[i]),y[i],sqrt(x[i]),z[i],col="blue")
    }
  }
  lines(sqrt(x),z,col="black",type="l")
  lines(sqrt(x),y,col="black",type="l")
  
  #x11()
  h = hist(distd[,(betti+1)],breaks = 50,plot = F)
  ccat = cut(h$breaks,c(-Inf,quantile(distd[,(betti+1)],0.95),Inf))
  plot(h,col=c("white","red")[ccat],xlab="red=pval. .blue=d-statis",
       main=mainhist)
  abline(v = quantile(distd[,(betti+1)],0.95),col="red")#p-val
  abline(v = d,col="blue")#statistic
}
#rechaza = reject in spanish
rechaza <- function(DMD,meanbc,distd,betti){
  #DMD: matrix with topological information
  #meanbc:  output of meanbc function
  #distd: distribution of distances to meanbc
  #betti: degree 0 or 1
  x = meanbc[,1]               #x axos
  y = meanbc[,(betti+2)]           #meanbck
  z = numeric(length(x))      #betti curve (not MBC)
  d = 0
  DMD1 = DMD[which(DMD[,1]==betti),2:3] #topological information
  for(i in 1:length(x)){
    z[i] = betti.fun(DMD1,x[i]) 
    d = d + abs(y[i]-z[i])
  }
  if(d>quantile(distd[,(betti+1)],0.95)){ #testing
    reject = 1
  }else{
    reject = 0
  }
  return(reject)
}

