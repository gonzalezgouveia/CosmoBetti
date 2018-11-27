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
#inicio de tiempos
t1 <- proc.time()
print('empieza proceso')

#set working directory
setwd("C:/Users/fafa/Desktop/proyecto_cosmologia/datos2/")

#leer los datos
NombresDatos <- list.files(pattern = "*.dat")
NombresVariables <- c()
for(nombre in NombresDatos){
  aux1 <- strsplit(nombre,"[.]")[[1]][1] #separa en punto
  aux2 <- strsplit(aux1,"_")[[1]][c(2,4)] #separa en underscore
  NombreVariable <- Reduce(function(...){paste(...,sep="")},aux2)
  print(NombreVariable)
  NombresVariables <- c(NombresVariables,NombreVariable)
  assign(NombreVariable,as.matrix(read.table(nombre)))
}

print('### datos leidos ###')
t2 <- proc.time()
print(t2-t1)

#calcular homología persistente
NombresDiag <- c()
for(nombre in NombresVariables){
  print(paste('calculando diagrama de: ',nombre))
  NombreDiag <- paste('Diag',nombre,sep="_")
  diag <- alphaShapeDiag(X = get(nombre), printProgress = TRUE)
  assign(NombreDiag,diag)
  NombresDiag <- c(NombresDiag,NombreDiag)
  print(paste('diagrama creado como: ',NombreDiag))
}
print('### diagramas de persistencia calculados ###')
t3 <- proc.time()
print(t3-t2)

#calcular curvas de Betti
NombresCurvas <- c()
for(nombre in NombresDiag){
  print(paste('calculando curva de betti de: ',nombre))
  auxnombre <- strsplit(nombre,"_")[[1]][2]
  NombreCurva <- paste('BC',auxnombre,sep="_")
  bc <- betticurves(xmax = 40,l=500,ASDIAG=get(nombre))
  assign(NombreCurva,bc)
  NombresCurvas <- c(NombresCurvas,NombreCurva)
  print(paste('curva de betti creada como: ',NombreCurva))
}
print('### curvas de betti calculadas ###')
t4 <- proc.time()
print(t4-t3)

#cambiar directorio donde guardare las imagenes
setwd("C:/Users/fafa/Desktop/proyecto_cosmologia/graficasreporte2")

#graficar de las curvas
# xmax = 40
# xseq <-seq(0,xmax,length.out = 500) 
# xlim_vec = c(10,25,40)
# for(k in 1:3){
#   print(paste('Curva de Betti de orden: ',k))
#   xlim = xlim_vec[k]
#   pdf(paste('todojunto_betti',k,'.pdf',sep=""),width = 6,height = 6)
#   grafica1(BC_F4Box1,k,'todos los modelos',add=FALSE,col = 1,xseq,xlim)
#   ColorIndex <- c(rep(1,4),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5))
#   i <- 1
#   for(nombre in NombresCurvas[-1]){
#     col <- ColorIndex[i]
#     i <- i+1
#     grafica1(get(nombre),k,'',add=TRUE, col = col,xseq,xlim)
#   }
#   legend("topright", inset=.05, title="Modelos",
#          c("F4","F5","F6","GR","N1","N5"), fill=c(1:6), horiz=F)
#   dev.off()
# }

#CALCULAR LAS CURVAS PROMEDIO PARA LOS 6 MODELOS
# NombresCurvasPromedio <- c()
# for(nombre in c('F4','F5','F6','GR','N1','N5')){
#   print(paste('trabajando con el modelo: ',nombre))
#   NombreLista <- paste('MBC',nombre,sep = '_')
#   assign(NombreLista,list())
#   NombresCurvasPromedio <- c(NombresCurvasPromedio,NombreLista)
#   for(k in 0:2){
#     print(paste('calculando MBC para grado de homologia: ',k))
#     NombreMBC <- paste('MBC_',k,'_',nombre,sep='')
#     assign(NombreMBC,rep(0,500))
#     for(index in 1:5){
#       NombreVar <- paste('BC_',nombre,'Box',index,sep="")
#       assign(NombreMBC,get(NombreMBC)+get(NombreVar)[[k+1]])
#     }
#     assign(NombreLista,c(get(NombreLista),list(get(NombreMBC)/5)))
#   }
# }
# print('Curvas Promedio calculadas')

#HACER LA GRAFICA CON LAS 6 LINEAS DE LAS PROMEDIOS
# xmax = 40
# xseq <-seq(0,xmax,length.out = 500) 
# xlim_vec = c(10,25,40)
# for(k in 1:3){
#   print(paste('Curva de Betti PROMEDIO de orden: ',k))
#   xlim = xlim_vec[k]
#   pdf(paste('promedio_betti',k,'.pdf',sep=""),width = 6,height = 6)
#   grafica1(MBC_F4,k,'todos los modelos',add=FALSE,col = 1,xseq,xlim)
#   #ColorIndex <- c(rep(1,4),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5))
#   i <- 1
#   for(nombre in NombresCurvasPromedio[-1]){
#     #col <- ColorIndex[i]
#     i <- i+1
#     grafica1(get(nombre),k,'',add=TRUE, col = i,xseq,xlim)
#   }
#   legend("topright", inset=.05, title="Modelos",
#          c("F4","F5","F6","GR","N1","N5"), fill=c(1:6), horiz=F)
#   dev.off()
# }
# print('grafica de los promedios')

#IMPRIMIR TODAS LAS DIFERENCIAS DOS A DOS DE LAS PROMEDIOS
# xmax = 40
# xseq <-seq(0,xmax,length.out = 500) 
# xlim_vec = c(10,25,40)
# 
# ModelosIndex <- c('F4','F5','F6','GR','N1','N5')
# for(mod1 in 1:5){
#   #asignar curva de betti promedio para modelo 1
#   NombreMBC1 <- paste('MBC_',ModelosIndex[mod1],sep="")
#   MBC1 <- get(NombreMBC1)
#   for(mod2 in (mod1+1):6){
#     #asignar curva de betti promedio para modelo 2
#     NombreMBC2 <- paste('MBC_',ModelosIndex[mod2],sep="")
#     MBC2 <- get(NombreMBC2)
#     for(k in 1:3){
#       #generando pdf
#       xlim <- xlim_vec[k]
#       m1 <- ModelosIndex[mod1];    m2 <- ModelosIndex[mod2]
#       titulo <- paste('diff_betti_mod',m1,'-',m2,'_betti_',k,sep="")
#       pdf(paste(titulo,'.pdf',sep=""),width = 6,height = 6)
#       #png(paste(titulo,'.png',sep=""),width = 480,height = 480)
#       bettiplotdiff(xseq,MBC1,MBC2,k,main=titulo
#                     ,add=F,col=1,xlim)
#       dev.off()
#     }
#   }
# }
# print('graficas graficadas :p')

#CALCULAR DESVIACION ESTANDAR DE LAS CURVAS
ModelosIndex <- c('F4','F5','F6','GR','N1','N5')
for(mod in 1:6){
  #print(paste('Trabajando con modelo: ',ModelosIndex[mod]))
  #inicializando variable para la desviacion estandar
  NombreVarVar <- paste('SD_MBC_',ModelosIndex[mod],sep="")
  assign(NombreVarVar,list())
  #print(paste('Creando Lista',NombreVarVar))
  #asignando MCB_"modelos" correspondiente
  NombreMBC <- paste('MBC_',ModelosIndex[mod],sep="")
  MCB <- get(NombreMBC)
  #recuerda que MCB_"modelos" es una lista de 3 vectores
  for(k in  1:3){
    #print(paste('Grafo de homologia: ',k,'. Para:',NombreVarVar))
    MBC_k <- MCB[[k]] #asignando la MCB corresp al grado de homologia
    var <- rep(0,500)
    for(l in 1:500){
      for(muestra in 1:5){
        NombreVarMuestra <- paste('BC_',ModelosIndex[mod],'Box',muestra,sep="")
        Lista <- get(NombreVarMuestra)
        Vector <- Lista[[k]]
        Valor <- Vector[l]
        aux <- 1/5*(Valor-MBC_k[l])^2 #calculo de varianza
        var[l] <- var[l]+aux
      }
      var[l] <- sqrt(var[l]) #esto es para tomar la desviacion
    }
    assign(NombreVarVar,c(get(NombreVarVar),list(var)))
    print(paste('calculada varianza para curva:',NombreVarVar))
  }
}
print('Calculadas todas las varianzas')

#IMPRIMIR DIFERENCIAS DOS A DOS CON DESVIACION ESTANDAR
xmax = 40
xseq <-seq(0,xmax,length.out = 500) 
xlim_vec = c(10,25,40)

ModelosIndex <- c('F4','F5','F6','GR','N1','N5')
for(mod1 in 1:5){
  #asignar curva de betti promedio para modelo 1
  NombreMBC1 <- paste('MBC_',ModelosIndex[mod1],sep="")
  MBC1 <- get(NombreMBC1)
  NombreSD1 <- paste('SD_MBC_',ModelosIndex[mod1],sep="")
  SD1 <- get(NombreSD1)
  for(mod2 in (mod1+1):6){
    #asignar curva de betti promedio para modelo 2
    NombreMBC2 <- paste('MBC_',ModelosIndex[mod2],sep="")
    MBC2 <- get(NombreMBC2)
    NombreSD2 <- paste('SD_MBC_',ModelosIndex[mod2],sep="")
    SD2 <- get(NombreSD2)
    for(k in 1:3){
      #generando pdf
      xlim <- xlim_vec[k]
      m1 <- ModelosIndex[mod1];    m2 <- ModelosIndex[mod2]
      titulo <- paste('diff_banda_mod',m1,'-',m2,'_betti_',k,sep="")
      #pdf(paste(titulo,'.pdf',sep=""),width = 6,height = 6)
      png(paste(titulo,'.png',sep=""),width = 480,height = 480)
      bettiplotbanda(xseq,MBC1,MBC2,SD1,SD2,k,main=titulo,xlim)
      dev.off()
    }
  }
}
print('graficas graficadas :p')
t5 <- proc.time()
print(t5-t4)

#PARA UNA SEGUNDA ENTREGA PROPONER SI HACER LA PRUEBA DE HIPOTESIS
#ESTO SERIA UNA TABLA CON LAS COMPARACIONES DOS A DOS
#DE ESTA TABLA REGRESARIA UN "RECHAZA QUE SON IGUALES" O NO

#tiempos
print('total')
print(t5-t1)


