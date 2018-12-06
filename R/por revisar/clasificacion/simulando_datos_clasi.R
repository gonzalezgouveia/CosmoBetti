######################################################################
#
#simulando los datos
#
######################################################################
n <- 50
B <- 20
#necesito crear los datos
for(i in 1:(B)){
  x <- sphereUnif(n,d=2)
  filename <- paste(data_directory,"sphere",i,sep="")
  output_file <- file(filename,'w')
  write(t(x),file=output_file,ncolumns=3,append=FALSE)
  close(output_file)
}
print("1 de 3")
for(i in 1:(B)){
  x <- rcubo(n,r=(pi/6)^(1/3))
  #x <- rcubo(n,r=0.9)
  filename <- paste(data_directory,"cube",i,sep="")
  output_file <- file(filename,'w')
  write(t(x),file=output_file,ncolumns=3,append=FALSE)
  close(output_file)
}
print("2 de 3")
for(i in 1:(B)){
  x <- klein.unif(k=n,sigma=0.0)
  for(j in 1:3){
    x[,j] <- x[,j]/(max(x[,j])-min(x[,j]))*2
  }
  filename <- paste(data_directory,"klein",i,sep="")
  output_file <- file(filename,'w')
  write(t(x),file=output_file,ncolumns=3,append=FALSE)
  close(output_file)
}
print("listo")