#Programa de analisis de datos de Cosmologia

#Se pretende analizar datos de galaxias de diferentes modelos a 
#traves de las curvas de Betti.
#Se calculan tres cosas:
  #1 Curvas de Betti para detectar diferencias
  #2 Curva de Betti promedio y error estandar
  #3 tiempos que tarde al algoritmo

#librerias
library(TDA) #libreria/modulo de R para análisis topológico de datos
library(ggplot2) #libreria/modulo de gráficas bonitas

#funciones
betti.fun <- function(bcmat,t){
  #return the number of birth - death
  #bcmat es una matriz
  nb <- sum(bcmat[,1]<t) #number of birth
  nd <- sum(bcmat[,2]<t) #number of death
  return(nb-nd)
}
betticurves <- function(xmax,l=500,ASDIAG){
  #esta funcion regresa una lista con las primeras tres
  #curvas de betti
  #xmax es el límite superior de las curvas de betti
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
bettilab <- function(k){
  if(k==1){
    ylab='número de componentes conexas'
  }else if(k==2){
    ylab='número de ciclos 2D'
  }else if (k==3){
    ylab='número de agujeros 3D'
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
  #k es el grado de homología puede ser 1,2,3
  bcdiff <- bc1[[k]][1,]-bc2[[k]][1,]
  if(add==F){
    plot(xseq,bcdiff,main=main,type='l',col=col,xlim=c(0,xlim),
         ylab=bettilab(k))
  }else{
    lines(xseq,bcdiff,main=main,type='l',col=col,xlim=c(0,xlim),
          ylab=bettilab(k))
  }
}
#inicio de tiempos
t1 <- proc.time()
print('empieza proceso')

#set working directory
setwd("C:/Users/fafa/Desktop/proyecto_cosmologia/datos2/")

#leer los datos
NombresDatos <- list.files(pattern = "*.dat")
for(nombre in NombresDatos){
  aux1 <- strsplit(nombre,"[.]")[[1]][1] #separa en punto
  aux2 <- strsplit(aux1,"_")[[1]][c(2,4)] #separa en underscore
  NombreVariable <- Reduce(function(...){paste(...,sep="")},aux2)
  print(NombreVariable)
  assign(NombreVariable,as.matrix(read.table(nombre)))
}
print('datos leidos')
t2 <- proc.time()
print(t2-t1)

#calcular homología persistente
Diag_gr1 <- alphaShapeDiag(X = grbox1, printProgress = TRUE)
Diag_gr2 <- alphaShapeDiag(X = grbox2, printProgress = TRUE)
Diag_n1 <- alphaShapeDiag(X = n1box1, printProgress = TRUE)
Diag_n5 <- alphaShapeDiag(X = n5box1, printProgress = TRUE)
print('diagramas de persistencia calculados')
t3 <- proc.time()
print(t3-t2)
#calcular curvas de Betti
bc_gr1 <- betticurves(xmax = 40,l=500,ASDIAG=Diag_gr1)
bc_gr2 <- betticurves(xmax = 40,l=500,ASDIAG=Diag_gr2)
bc_n1 <- betticurves(xmax = 40,l=500,ASDIAG=Diag_n1)
bc_n5 <- betticurves(xmax = 40,l=500,ASDIAG=Diag_n5)
print('curvas de betti calculadas')
t4 <- proc.time()
print(t4-t3)

#graficar de las curvas
xmax = 40
xseq <-seq(0,xmax,length.out = 500) 
xlim_vec = c(10,25,40)
for(k in 1:3){
  xlim = xlim_vec[k]
  pdf(paste('todojunto_betti',k,'.pdf',sep=""),
      width = 6,height = 6)
  grafica1(bc_gr1,k,'todos los modelos',add=FALSE,col = 1,xseq,xlim)
  grafica1(bc_gr2,k,'',add=TRUE, col = 4,xseq,xlim)
  grafica1(bc_n1,k,'',add = TRUE,col=2,xseq,xlim)
  grafica1(bc_n5,k,'',add=T,col=3,xseq,xlim)
  legend("topright", inset=.05, title="Modelos",
         c("gr1","gr2","n5","n1"), fill=c(1,4,2,3), horiz=F)
  dev.off()
}


for(k in 1:3){
  xlim = xlim_vec[k]
  pdf(paste('diferencias_betti',k,'.pdf',sep=""),
      width = 6,height = 6)
  bettiplotdiff(xseq,bc_gr1,bc_n1,k,main=paste('Diferencias curvas de Betti. beta = ',(k-1))
                ,add=F,col=1,xlim)
  bettiplotdiff(xseq,bc_gr1,bc_gr2,k,main=paste('')
                ,add=T,col=2,xlim)
  legend("bottomright", inset=.05, title="Modelos",
         c("gr1-gr2","gr1-n1"), fill=c(2,1), horiz=F)
  dev.off()
}
print('graficas graficadas :p')
t5 <- proc.time()
print(t5-t4)
#hacer graficas con diferencias
  #una idea puede ser hacer gr1 como el estandar
  #y las demas 3 gráficas serian las diferencias con respecto a gr1
  #de esta forma, se vería cómo gr2 y n5 son cercanas a gr1
  #y como n1 es más lejana a gr1...
  #la próxima misión sería determinar si esta diferencia es
  #significativa, ¿qué tamaño de muestra es necesario para esto?
  #la cuestions es tener n muestras de gr1 y n muestras de n1
  #luego hacer tomar la estadistica de las diferencias y determinar si
  #son diferentes
  #propongo un experimento piloto con 10 de cada una
  #se expondría una representación gráfica
    #con curvas de betti promedio y bandas

#tiempos
print('total')
print(t5-t1)


