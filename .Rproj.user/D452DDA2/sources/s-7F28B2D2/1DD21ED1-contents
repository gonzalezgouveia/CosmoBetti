# funciones que tienen que ver con las graficas

library(TDA) #libreria/modulo de R para an?lisis topol?gico de datos
library(ggplot2) #libreria/modulo de gr?ficas bonitas

bettilab <- function(k){
  if(k==1){
    ylab='n?mero de componentes conexas'
  }else if(k==2){
    ylab='n?mero de ciclos 2D'
  }else if (k==3){
    ylab='n?mero de agujeros 3D'
  }
  return(ylab)
}
grafica1 <- function(ph_list,k,main,add=FALSE,col,xseq,xlim,ylab){
  if(add==FALSE){
    plot(xseq,ph_list[[k]][1,],type='S',main = paste('betti = ',(k-1),main),
         col=col,xlim=c(0,xlim),ylab=bettilab(k))
  }else{
    lines(xseq,ph_list[[k]][1,],type='S',main = paste('betti = ',(k-1),main),
          col=col,xlim=c(0,xlim),ylab=bettilab(k))
  }
}
bettiplotdiff <- function(xseq,bc1,bc2,k,main='',add,col,xlim,ylab){
  #k es el grado de homolog?a puede ser 1,2,3
  bcdiff <- bc1[[k]][1,]-bc2[[k]][1,]
  if(add==F){
    plot(xseq,bcdiff,main=main,type='l',col=col,xlim=c(0,xlim),
         ylab=bettilab(k))
  }else{
    lines(xseq,bcdiff,main=main,type='l',col=col,xlim=c(0,xlim),
          ylab=bettilab(k))
  }
}

bettiplotbanda <- function(xseq,bc1,bc2,sd1,sd2,k,main='',xlim,ylab){
  bcdiff <- bc1[[k]][1,]-bc2[[k]][1,]
  #curva para modelo 2
  plot(xseq,bcdiff,main=main,type='l',col=1,xlim=c(0,xlim),
       ylab=bettilab(k))
  #desviacion para modelo 2
  sd2sup <- bcdiff+sd2[[k]]
  sd2inf <- bcdiff-sd2[[k]]
  lines(xseq,sd2sup,col=2)
  lines(xseq,sd2inf,col=2)
  polygon(c(xseq, rev(xseq)), c(sd2sup, rev(sd2inf)),
          col = rgb(1, 0, 0,0.1), border = NA)
  lines(xseq,bcdiff,col=1,lwd=2)
  #curva para modelo 1
  abline(h=0,col=1)  #eje x
  #banda para modelo 1
  sd1sup <- sd1[[k]]
  sd1inf <- -sd1[[k]]
  lines(xseq,sd1sup,col=4)
  lines(xseq,sd1inf,col=4)
  polygon(c(xseq, rev(xseq)), c(sd1sup, rev(sd1inf)),
          col = rgb(0, 0, 1,0.1), border = NA)
  abline(h=0,col=1,lwd=2)
}
