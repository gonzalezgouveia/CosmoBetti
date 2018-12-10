# funciones que tienen que ver con el vector de betti

library(TDA) #libreria/modulo de R para an?lisis topol?gico de datos
library(tidyverse)

#funciones
##########

#' compute the difference between birth-death at time t from a matrix.
#' internal function used in bettimatrix
#'
#' @param bcmat is a ?times? matrix
#' @param t is the time or distance to compute betti.fun
#' @return the number of birth - death
betti.fun <- function(bcmat,t){
  nb <- sum(bcmat[,1]<t)
  nd <- sum(bcmat[,2]<t)
  return(nb-nd)
}

#' Compute Betti Curves
#'
#' Compute Betti Curves until some order from persistence diagram
#' if order=-1 then the order is infered from diag
#' if xmax=-1 the max lenght is infered from diag
#'
#' @param diag is a persistence diagram from TDA o TDAstats
#' @param order is the maximum order the betti curve is computed
#' @param xmax is the maximum length/time the betti curve is computed
#' @param l is the dimension of each betticurve
#' @return a 3 times l matrix with the betticurve
#' @examples
#' library(TDA)
#' data <- cbind(runif(10), runif(10))
#' alpha_complex <- alphaComplexDiag(data)
#' alpha_diag <- alpha_complex$diagram
#' bettimatrix(alpha_complex)
bettimatrix <- function(diag, order = -1, xmax=-1, l = 100){
  # check if diag is a matrix
  # PENDIENTE: HACER FUNCION PARA CHECAR QUE ES UN DIAG VALIDO
  if (is.matrix(diag)) {
  } else if (is.list(diag)) {
    diag <- diag$diagram
  } else {
    print('diag not a valid persistence matrix')
    break
  }

  # assign xmax automatically, info in description
  if (xmax == -1){
    xmax <- max(diag[-1,2:3]) * 1.05 # maximum plus 5%
  }

  # assign order automatically, info in description
  if (order == -1){
    order <- max(diag[,1]) # starts at 0
  }

  # create the x dimension for the bettimatrix
  xseq <- seq((xmax/(2*l)), xmax, length.out = l)

  bettimatrix <- matrix(NA, ncol = 3, nrow = 0)
  colnames(bettimatrix) <- c('order', 'xseq', 'value')

  for(i in 0:order){
    if(order==0){
      phx <- diag[which(diag[,1]==i), 2:3][-1,]
    } else {
      phx <- diag[which(diag[,1]==i), 2:3]
    }
    for(j in 1:l){
      bettivalue <- betti.fun(phx, xseq[j])
      bettimatrix <- rbind(bettimatrix,
                           c(i, xseq[j], bettivalue))
    }
  }
  bettimatrix <- as_tibble(bettimatrix)
  bettimatrix$order <- factor(bettimatrix$order)
  return(bettimatrix)
}
