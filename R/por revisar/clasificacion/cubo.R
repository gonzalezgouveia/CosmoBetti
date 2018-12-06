#simular de un cubo
n=10
rcubo <- function(n,r=1,d=3){
  nube <- numeric(0)
  for(i in 1:n){
    pos <- sample(d,1)
    signo <- sample(c(-r,r),1)
    punto <- numeric(d)
    punto[pos] <- signo
    punto[-pos] <- runif((d-1),min=-r,max=r) 
    nube <- rbind(nube,punto)
  }
  return(nube)
}
#rcubo(50,r=(pi/6)^(1/3))
