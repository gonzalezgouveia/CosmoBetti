#hacer algoritmo de clasificacion con curvas de betti

######################################################################
#
#parametros iniciales
#
######################################################################

#asignando directorio de la muestra
#sample_directory <- choose.dir()
sample_directory <- "C:\\Users\\fafa\\Desktop\\bubenik\\TDA_bubenik_tesis"
setwd(sample_directory) # sample directory for windows users

input <- "PC" # PC or Mac - for linux use Mac
data_directory <- "nonrigid3d/"
#perseus_directory <- "perseus_output/"
pl_directory <- "pl_output/"

grid_size <- 1 # grid used for Perseus and Persistence Landscape Toolbox
max_depth <- 200 # max number of persistence landscape functions
max_filtration_value <- grid_size
pl_param_vals <- seq(0,grid_size,0.01) # values at which to evaluate the persistence landscape

max_pose_num <- 20 # 10 # maximum pose number for one figure, starts at 0
figures_to_use <- c("cube","klein")
figures_to_compare <- c("cube","klein")

N <- 100 # number of repeats in permutation test
cost <- 10 # lower: softer margin, higher: harder margin, default=1
num_folds <- 10 # number of folds for cross validation
show_individual_figures <- F # draw filtered simplicial complexes, persistence diagrams, persistence landscapes
pca_coords_to_use <- 1:10 # number of pca coordinates to use for classification using pca
degrees_to_use <- 0:2 # which homological degrees are used for classification


library(rgl) # needed for 3d graphicss
library(FNN) # needed for k-nearest neighbor
library(scatterplot3d)
library(e1071) # needed for PCA, SVM
library(kernlab)
library(TDA)

source("tda_functions.R") # load Peter Bubenik's topological data analysis functions.

######################################################################
#
#simulando los datos
#
######################################################################
#n <- 100
#B <- 500
##necesito crear los datos
#for(i in 1:(B)){
#  x <- sphereUnif(n,d=2)
#  filename <- paste(data_directory,figures_to_use[1],i,sep="")
#  output_file <- file(filename,'w')
#  write(t(x),file=output_file,ncolumns=3,append=FALSE)
#  close(output_file)
#}
#print("1 de 3")
#for(i in 1:(B)){
#  x <- rcubo(n,r=(pi/6)^(1/3))
#  filename <- paste(data_directory,figures_to_use[2],i,sep="")
#  output_file <- file(filename,'w')
#  write(t(x),file=output_file,ncolumns=3,append=FALSE)
# close(output_file)
#}
#print("2 de 3")
#for(i in 1:(B)){
#  x <- klein.unif(k=n,sigma=0.0)
#  for(j in 1:3){
#    x[,j] <- x[,j]/(max(x[,j])-min(x[,j]))*2
#  }
#  filename <- paste(data_directory,figures_to_use[3],i,sep="")
#  output_file <- file(filename,'w')
#  write(t(x),file=output_file,ncolumns=3,append=FALSE)
#  close(output_file)
#}
#
######################################################################
#
#calcular homologia, curvas de betti, vectores de curvas de betti
#
######################################################################
betti.fun <- function(bcmat,t){
  #return the number of birth - death
  #bcmat es una matriz
  nb <- sum(bcmat[,1]<t) #number of birth
  nd <- sum(bcmat[,2]<t) #number of death
  return(nb-nd)
}
# get list of figures
figures <- vector() # initialize vector
for (figure_name in figures_to_use) # for each type of figure
  for (i in 1:max_pose_num)
    if (file.exists(paste(data_directory,figure_name,i,sep="")))
      figures <- rbind(figures, c(figure_name,i))
num_figures <- dim(figures)[1]

l <-  500 # tamano del vector de curvas de betti
plm0 <- matrix(0,ncol = l,nrow = num_figures)
plm1 <- matrix(0,ncol = l,nrow = num_figures)
plm2 <- matrix(0,ncol = l,nrow = num_figures)

pose <- 1 # pose number between 1 and num_figures (148)
for (pose in 1:num_figures){ # for each pose
  # read the data
  file <- paste(figures[pose,1],figures[pose,2],sep="")
  
  print(file)
  x <- read.table(paste(data_directory,file,sep=""))
  diag <- alphaShapeDiag(x)$diag
  ph0 <- diag[which(diag[,1]==0),2:3][-1,]
  ph1 <- diag[which(diag[,1]==1),2:3]
  ph2 <- diag[which(diag[,1]==2),2:3]
  if(length(ph2)==2){ #cuando hay una sola componente en H2
    ph2 <- t(as.matrix(ph2))
  }
  
  #aqui quiero hacer la cruva de betti para ph0 y para ph1
  xseq <- seq(0,1.1,length.out = l)
  betticurve0 <- numeric(l)
  betticurve1 <- numeric(l)
  betticurve2 <- numeric(l)
  for(j in 1:l){
    betticurve0[j] <- betti.fun(ph0,xseq[j])
    betticurve1[j] <- betti.fun(ph1,xseq[j])
    betticurve2[j] <- betti.fun(ph2,xseq[j])
  }
  plm0[pose,] <- betticurve0
  plm1[pose,] <- betticurve1
  plm2[pose,] <- betticurve2
}
plm_list <- list(plm0,plm1,plm2)
#graficando las curvas de betti promedio para las figuras
#################################################################
#agregar BANDAS EN ESTA SECCION
#observar que plm0 es el vector con todas las curvas de betti
#sacar las desviacion de estas matrices
#################################################################
library(ggplot2)
B = num_figures/2
#construyendo MBC
mbc0figA <- colMeans(plm0[1:B,])
mbc0figB <- colMeans(plm0[(B+1):(2*B),])
mbc1figA <- colMeans(plm1[1:B,])
mbc1figB <- colMeans(plm1[(B+1):(2*B),])
mbc2figA <- colMeans(plm2[1:B,])
mbc2figB <- colMeans(plm2[(B+1):(2*B),])

upp0figA <- numeric()
down0figA <- numeric()
upp0figB <- numeric()
down0figB <- numeric()
upp1figA <- numeric()
down1figA <- numeric()
upp1figB <- numeric()
down1figB <- numeric()
upp2figA <- numeric()
down2figA <- numeric()
upp2figB <- numeric()
down2figB <- numeric()
#construyendo bandas de confianza
for(j in 1:l){
  #betti0
  upp0figA[j] <- mbc0figA[j]+sd(plm0[1:B,j])
  down0figA[j] <- mbc0figA[j]-sd(plm0[1:B,j])
  upp0figB[j] <- mbc0figB[j]+sd(plm0[(B+1):(2*B),j])
  down0figB[j] <- mbc0figB[j]-sd(plm0[(B+1):(2*B),j])
  #betti1
  upp1figA[j] <- mbc1figA[j]+sd(plm1[1:B,j])
  down1figA[j] <- mbc1figA[j]-sd(plm1[1:B,j])
  upp1figB[j] <- mbc1figB[j]+sd(plm1[(B+1):(2*B),j])
  down1figB[j] <- mbc1figB[j]-sd(plm1[(B+1):(2*B),j])
  #betti2
  upp2figA[j] <- mbc2figA[j]+sd(plm2[1:B,j])
  down2figA[j] <- mbc2figA[j]-sd(plm2[1:B,j])
  upp2figB[j] <- mbc2figB[j]+sd(plm2[(B+1):(2*B),j])
  down2figB[j] <- mbc2figB[j]-sd(plm2[(B+1):(2*B),j])
}
data1 = as.data.frame(cbind(xseq,mbc0figA,upp0figA,down0figA,
                            mbc0figB,upp0figB,down0figB,
                            mbc1figA,upp1figA,down1figA,
                            mbc1figB,upp1figB,down1figB,
                            mbc2figA,upp2figA,down2figA,
                            mbc2figB,upp2figB,down2figB)) 

#para betti0
p1 <- ggplot(data1,aes(x=xseq))+
  #ejes y titulo
  xlab("distancia")+ylab("cantidad")+
  labs(title=paste("A=",figures_to_use[1],"B=",figures_to_use[2]))+
  xlim(c(0,0.5))+
  #lineas
  geom_line(aes(y=mbc0figA,colour="A"))+
  geom_line(aes(y=mbc0figB,color="B"))+
  #banda para figura A
  geom_ribbon(aes(x=xseq,ymax = upp0figA,ymin= down0figA),fill="green",alpha=0.2,show.legend = F)+
  geom_line(aes(y=upp0figA,colour = "A"))+
  geom_line(aes(y=down0figA,colour = "A"))+
  #banda para figura B
  geom_ribbon(aes(x=xseq,ymax = upp0figB,ymin= down0figB),fill="red",alpha=0.2,show.legend = F)+
  geom_line(aes(y=upp0figB,colour = "B"))+
  geom_line(aes(y=down0figB,colour = "B"))+
  #escogiendo colores de graficas
  scale_color_manual("",breaks=c("A","B"),values=c("A"="green","B"="red"))
print(p1)
#para betti1
p2 <- ggplot(data1,aes(x=xseq))+
  #ejes y titulo
  xlab("distancia")+ylab("cantidad")+
  labs(title=paste("A=",figures_to_use[1],"B=",figures_to_use[2]))+
  xlim(c(0,1))+
  #lineas
  geom_line(aes(y=mbc1figA,colour="A"))+
  geom_line(aes(y=mbc1figB,color="B"))+
  #banda para figura A
  geom_ribbon(aes(x=xseq,ymax = upp1figA,ymin= down1figA),fill="green",alpha=0.2,show.legend = F)+
  geom_line(aes(y=upp1figA,colour = "A"))+
  geom_line(aes(y=down1figA,colour = "A"))+
  #banda para figura B
  geom_ribbon(aes(x=xseq,ymax = upp1figB,ymin= down1figB),fill="red",alpha=0.2,show.legend = F)+
  geom_line(aes(y=upp1figB,colour = "B"))+
  geom_line(aes(y=down1figB,colour = "B"))+
  #escogiendo colores de graficas
  scale_color_manual("",breaks=c("A","B"),values=c("A"="green","B"="red"))
print(p2)
#para betti2
p3 <- ggplot(data1,aes(x=xseq))+
  #ejes y titulo
  xlab("distancia")+ylab("cantidad")+
  labs(title=paste("A=",figures_to_use[1],"B=",figures_to_use[2]))+
  xlim(c(0,1.1))+
  #lineas
  geom_line(aes(y=mbc2figA,colour="A"))+
  geom_line(aes(y=mbc2figB,color="B"))+
  #banda para figura A
  geom_ribbon(aes(x=xseq,ymax = upp2figA,ymin= down2figA),fill="green",alpha=0.2,show.legend = F)+
  geom_line(aes(y=upp2figA,colour = "A"))+
  geom_line(aes(y=down2figA,colour = "A"))+
  #banda para figura B
  geom_ribbon(aes(x=xseq,ymax = upp2figB,ymin= down2figB),fill="red",alpha=0.2,show.legend = F)+
  geom_line(aes(y=upp2figB,colour = "B"))+
  geom_line(aes(y=down2figB,colour = "B"))+
  #escogiendo colores de graficas
  scale_color_manual("",breaks=c("A","B"),values=c("A"="green","B"="red"))
print(p3)

#x11()
#B = num_figures/2
#plot(xseq,colMeans(plm0[1:B,]),type="l",main=paste("h0",figures_to_use[1]))
#plot(xseq,colMeans(plm0[(B+1):(2*B),]),type="l",main=paste("h0",figures_to_use[2]))
#plot(xseq,colMeans(plm1[1:B,]),type="l",main=paste("h1",figures_to_use[1]))
#plot(xseq,colMeans(plm1[(B+1):(2*B),]),type="l",main=paste("h1",figures_to_use[2]))
#plot(xseq,colMeans(plm2[1:B,]),type="l",main=paste("h2",figures_to_use[1]))
#plot(xseq,colMeans(plm2[(B+1):(2*B),]),type="l",main=paste("h2",figures_to_use[2]))

#rm(plm0,plm1)

#save to R Data file
#save(figures,num_figures,plm_list,grid_size,max_depth,pl_param_vals,file="nonrigid3d_landscape.RData")

######################################################################
#
#modelo de clasificacion con svm
#
######################################################################
#corriendo solo la parte de svm de bubenik
require("kernlab")
pl <- numeric()
for (i in degrees_to_use)
  pl <- cbind(pl, plm_list[[i+1]])
#pl <- pl[1:500,]
#pl <- pl[c(1:250,501:750),]
#pl <- pl[251:750,]
sample_classes <- figures[,1]

classes_to_use <- figures_to_use

indices_to_use <- c()
for (i in 1:length(classes_to_use))
  indices_to_use <- c(indices_to_use,which(sample_classes[]==classes_to_use[i]))

svm_classes <- factor(sample_classes[indices_to_use,drop=F])

print("Classification using SVM and k-fold cross validation")
classification_svm(pl,svm_classes,num_folds,cost)