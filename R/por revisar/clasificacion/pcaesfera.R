#haciendo pca a los datos de cubos y esferas
library(ElemStatLearn)
library(mvtnorm)
library(TDA)
n <- 50 #cantidad de puntos de la esfera
B <- 500 #cantidad de nubes de datos. Solo poner cantidades pares
dat1 <- numeric(0)
for(i in 1:(B/2)){
  #a <- sphereUnif(n,2,r=1)
  a <- rcubo(n,r=(pi/6)^(1/3))
  #a <- a[order(a[,1]),]
  dat1 <- rbind(dat1,matrix(a,nrow=1))
}
dat2 <- numeric(0)
for(i in 1:(B/2)){
  #a <- torusUnif(n,1,3)
  a <- klein.unif(k=n,sigma=0.0)
  for(j in 1:3){
    a[,j] <- a[,j]/(max(a[,j])-min(a[,j]))*2
  }
  #a <- rcubo(n,r=(pi/6)^(1/3))
  #a <- sphereUnif(n,2)
  #a <- a[order(a[,1]),]
  dat2 <- rbind(dat2,matrix(a,nrow=1))
}
dat <- rbind(dat1,dat2)
y <- c(rep(0,B/2),rep(1,B/2))

train <- dat
label <- as.factor(y)
covtrain <- cov(train)
#pca
train_pc <- prcomp(covtrain)
train_pc2 <- prcomp(cor(train))
train_pc <- prcomp(train)
varex <- train_pc$sdev^2/sum(train_pc$sdev^2)
varcum <- cumsum(varex)
result <- data.frame(num=1:length(train_pc$sdev),
                     ex=varex,
                     cum=varcum)

#plot(result$num,result$cum,type="b",xlim=c(0,150),
#     main="Varianza explicadad por las primeras 150 Componentes",
#     sub="caso esfera-cubo",
#     xlab="Numero de compoenentes",ylab="Varianza explicada")
#abline(v=3,lty=2)

train_score <- as.matrix(train) %*% train_pc$rotation[,1:25]
train <- cbind(label,as.data.frame(train_score))

colors <- rainbow(length(unique(train$label)))
names(colors) <- unique(train$label)
plot(train$PC1,train$PC2,type="n",main="Primeras dos componentes principales",
     sub="caso cubo-klein",xlab="primera componente",ylab="segunda componente")
text(train$PC1,train$PC2,label=train$label,col=colors[train$label])
legend("topleft",c("(0) cubo","(1) klein"),col=c(2,5),lty=1,lwd=3)
