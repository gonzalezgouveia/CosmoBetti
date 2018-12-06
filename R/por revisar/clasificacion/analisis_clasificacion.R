library(rgl)
n=50
#x <-  rcubo(n,r=(pi/6)^(1/3))
x <-  rcubo(n,r=.9)
y <- klein.unif(k=n,sigma=0.0)
for(j in 1:3){
  y[,j] <- y[,j]/(max(y[,j])-min(y[,j]))*2
}
z <- sphereUnif(n,d=2)
plot3d(x)
#plot3d(x,add=T,col=2)
diagx <- alphaShapeDiag(x)$diagram
diagy <- alphaShapeDiag(y)$diagram
diagz <- alphaShapeDiag(z)$diagram

plot(diagx,c(0,1),main="cube")
plot(diagy,c(0,1),main="klein")
plot(diagz,c(0,1),main="sphere")


#analizando el cubo
ph2x <- diagx[which(diagx[,1]==2),2:3]

#analizando la sphere
ph2 <- diagz[which(diagz[,1]==2),2:3]
ph2[1,1]==ph2[1,2]
ph2[2,1]-ph2[2,2]

betticurve2 <- numeric(l)
for(j in 1:l){
  betticurve2[j] <- betti.fun(ph2,xseq[j])
}
plot(xseq,betticurve2,type="l")
