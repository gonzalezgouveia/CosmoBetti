library(spatstat)
library(TDA)

#need to run topo_function.R first

#simulating data
squ <- 1; beta <- 150;gamma <- 0.5;R <- 0.05
X <- rStrauss(beta,gamma,R,square(squ))
#n <- X$n
mydata = matrix(c(X$x,X$y),ncol=2)

#binomial process
n <- 100 
mydata = matrix(runif(2*n)*squ,ncol=2)

#compute this dirst, takes 5 to 10 minuts. Only compute it once.
meanbc = meanbc.fun(xmax = 0.25,n=100,B=1000) #MBC
#distd = distd.fun(meanbc,B=5000) #distribtion of distances
#hist(distd[,2],breaks = 50)

#correr esta funciona para la imagen 2.6 de la tesis
#en mydata va el proceso al que se le quiere sacar la imagen
graficabetti <- function(mydata,meanbc,betti,xlim){
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
  if(betti == 1){
    main = "Loops. Betti 1"
    ylab = "mean number of loops"
    ylim = 12
  }else if(betti == 0){
    main = "Components. Betti 0"
    ylab = "mean number of connected components"
    ylim = 100
  }
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
  #h = hist(distd[,(betti+1)],breaks = 50,plot = F)
  #ccat = cut(h$breaks,c(-Inf,quantile(distd[,(betti+1)],0.95),Inf))
  #plot(h,col=c("white","red")[ccat],xlab="red=pval. .blue=d-statis",
  #     main=mainhist)
  #abline(v = quantile(distd[,(betti+1)],0.95),col="red")#p-val
  #abline(v = d,col="blue")#statistic
}

#plot test for one process, mydata
graficabetti(mydata,meanbc,betti = 0,xlim = 0.12)
graficabetti(mydata,meanbc,betti = 1,xlim = 0.18)
#x11()
