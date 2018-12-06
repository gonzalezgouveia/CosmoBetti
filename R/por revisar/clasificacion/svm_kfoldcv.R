#tratando de usar ksvm para clasificar mis imagenes

#leyendo datos de la carpeta /nongrid3d
dat <- numeric(0)
for (pose in 1:num_figures){ # for each pose
  # read the data
  file <- paste(figures[pose,1],figures[pose,2],sep="")
  #print(file)
  x <- read.table(paste(data_directory,file,sep=""))
  dat <- rbind(dat,matrix(as.matrix(x),nrow=1))
}
#hist(dat[1,])
#a1 <- matrix(dat[2,],ncol=3)
#plot3d(a1)
#plot(alphaShapeDiag(a1)$diag,c(0,1))
#a2 <- matrix(dat[251,],ncol=3)
#plot3d(a2)
#plot(alphaShapeDiag(a2)$diag,c(0,1))

M <- dat
require("kernlab") 
num_folds <- min(num_folds,length(svm_classes))
# Perform k fold cross validation
folds <- cut(seq(1,nrow(M)),breaks=num_folds,labels=FALSE) # Create k equally size folds
folds <- folds[sample(length(folds))]  # Randomly shuffle the folds
# SVM for PL matrix - use ksvm from kernlab (slow); svm from e1071 results in stack overflow 
num_sv <- 0 # total number of support vectors
pred_table <- matrix(0,nrow=nlevels(svm_classes),ncol=nlevels(svm_classes))
for(i in 1:num_folds){
  # Segement your data by fold using the which() function 
  testIndices <- which(folds==i,arr.ind=TRUE)
  testSet <- M[testIndices,,drop=FALSE] # prevent 1 row matrix from turning into a vector
  trainSet <- M[-testIndices,,drop=FALSE]
  # Fit the model
  svm_model <- ksvm(trainSet,svm_classes[-testIndices],type="C-svc",scaled=c(),kernel="vanilladot",C=cost)
  num_sv <- num_sv + nSV(svm_model)
  # Make the prediction (the dependent variable, without using the Type)
  svm_pred <- predict(svm_model, testSet)
  pred_table <- pred_table + table(pred = svm_pred, true = svm_classes[testIndices])
}
print(pred_table)
pred_accuracy <- sum(diag(pred_table)) / sum(pred_table)
print(paste("Average number of Support Vectors:",num_sv/num_folds))
print(paste("Predication Accuracy:",pred_accuracy)) 

#########
#tratando de usar kernels para clasificar mis datos
Gram <- crossprod(t(M))
classification_svm_gram(Gram,svm_classes,num_folds,cost)
