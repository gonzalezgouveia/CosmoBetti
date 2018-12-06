#resultados experimento de clasificacion
errores <- c(27.2,23.6,11.2)
metodos <- c("LDA","SVM","Bubenik")
plot(errores,xaxt="n",type="b",xlim=c(0.5,3.5),ylim=c(5,30),
     ylab="porcentaje",xlab="",main="Error de clasificacion")
axis(1,1:3,labels = F)
text(1:3+0.1, par("usr")[3] - 3, labels = as.vector(metodos), 
     srt = 50, pos = 2, xpd = T,cex=1)
