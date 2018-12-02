# funciones que tienen que ver con el vector de betti

library(TDA) #libreria/modulo de R para an?lisis topol?gico de datos
library(ggplot2) #libreria/modulo de gr?ficas bonitas

#funciones
betti.fun <- function(bcmat,t){
  #return the number of birth - death
  #bcmat es una matriz ¿de cuanto por cuanto?
  # que es la t
  nb <- sum(bcmat[,1]<t) #number of birth
  nd <- sum(bcmat[,2]<t) #number of death
  return(nb-nd)
}
betticurves <- function(xmax,l=500,ASDIAG){
  #esta funcion regresa una lista con las primeras tres
  #curvas de betti
  #xmax es el l?mite superior de las curvas de betti
  #(aqui se asumen los tres iguales)
  #l es la cantidad de elementos por vector
  #diag es el diagrama de persistencia
  diag <- ASDIAG$diag
  xseq <- seq(0,xmax,length.out = l)
  plm0 <- matrix(0,ncol = l,nrow = 1)
  plm1 <- matrix(0,ncol = l,nrow = 1)
  plm2 <- matrix(0,ncol = l,nrow = 1)
  betticurve0 <- numeric(l)
  betticurve1 <- numeric(l)
  betticurve2 <- numeric(l)
  ph0 <- diag[which(diag[,1]==0),2:3][-1,]
  ph1 <- diag[which(diag[,1]==1),2:3]
  ph2 <- diag[which(diag[,1]==2),2:3]
  for(j in 1:l){
    betticurve0[j] <- betti.fun(ph0,xseq[j])
    betticurve1[j] <- betti.fun(ph1,xseq[j])
    betticurve2[j] <- betti.fun(ph2,xseq[j])
  }
  plm0[1,] <- betticurve0
  plm1[1,] <- betticurve1
  plm2[1,] <- betticurve2

  plm_list <- list(plm0,plm1,plm2)
  return(plm_list)
}
