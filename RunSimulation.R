list.of.packages <- c("Matrix", "Rcpp", "RcppArmadillo","mvtnorm", "optimization", "matrixStats", "pracma", "penalizedLDA", "caret", "randomForest", "LiblineaR", "sparseLDA", "caret", "keras", "Rfast", "nloptr", "inline", "lbfgs")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(optimization)
library(matrixStats)
library(pracma)
#library(VaDA)
library(penalizedLDA)
library(caret)
library(randomForest)
library(LiblineaR)
library(sparseLDA)
library(caret)
library(keras)
sourceCpp("GPDAnonStatFinal.cpp")
source("GPDARversionFinal.R")
source("CNNcode.R")

n.varbs <- 5000
nTrain <- 100
nTest <- 500
delta <- 1
n.reps <- 50

#takes values from {1,2,3,4}
sim.set <- 1

if(sim.set == 1){
  
  alpha <- log(n.varbs)
  beta <- 3
  zeta.init <- 3
  
}

if(sim.set == 2){
  
  alpha <- log(n.varbs)
  beta <- 3
  zeta.init <- 3
}

if(sim.set == 3){
  
  alpha <- log(n.varbs)
  beta <- 3
  zeta.init <- 3
}

if(sim.set == 4){
  
  alpha <- log(n.varbs)
  beta <- 3
  zeta.init <- 0
  M <- matrix(0.95,nrow=5000,ncol=5000)
  diag(M) <- rep(1,5000)
  cholM <- t(chol(M))
}


Lambda_vec<-10^(seq(-3,2,1))
Lambda2_vec<-10^(seq(-3,2,1))
LambdaSparseLDA_vec<-c(1e-6,1e-4,1e-3,.1,1,10)
CNNsettings <- expand.grid(c(5,10,15,20), c(10,20,30,40))
names(CNNsettings) <- c("FILTERS","KERNELS")

Errors1 <- rep(0,n.reps)
Errors2 <- rep(0,n.reps)
Errors3 <- rep(0,n.reps)
Errors4 <- rep(0,n.reps)
Errors5 <- rep(0,n.reps)
Errors6 <- rep(0,n.reps)
Errors7 <- rep(0,n.reps)
Errors8 <- rep(0,n.reps)
Errors9 <- rep(0,n.reps)
Errors10 <- rep(0,n.reps)

RecordChoose1 <- matrix(0, n.reps, n.varbs)
RecordChoose2 <- matrix(0, n.reps, n.varbs)
RecordChoose3 <- matrix(0, n.reps, n.varbs)
RecordChoose4 <- matrix(0, n.reps, n.varbs)
RecordChoose5 <- matrix(0, n.reps, n.varbs)
RecordChoose6 <- matrix(0, n.reps, n.varbs)
RecordChoose8 <- matrix(0, n.reps, n.varbs)
RecordChoose9 <- matrix(0, n.reps, n.varbs)

MCC1 <- rep(0,n.reps)
MCC2 <- rep(0,n.reps)
MCC3 <- rep(0,n.reps)
MCC4 <- rep(0,n.reps)
MCC5 <- rep(0,n.reps)
MCC6 <- rep(0,n.reps)
MCC8 <- rep(0,n.reps)
MCC9 <- rep(0,n.reps)

for(sim.rep in 1:n.reps){
  
  source("GenSim.R")
  
  #Fit GPDA model
  GPobj <- GPDA.sparse.NonStat(mXtrain = X.train, vytrain = vy.train, mXtest = X.test, train.cycles = 5, test.cycles = 3, delta = 1, alpha = alpha, beta=beta)
  Errors1[sim.rep] <- sum(abs(round(c(GPobj$xi))!=vy.test))
  RecordChoose1[sim.rep,] <- round(GPobj$vw)
  MCC1[sim.rep] <- Calc.MCC(predvec =RecordChoose1[sim.rep,], truevec = gamma.true)
  
  #Fit VNPDA model(run if you installed VaDA)
  #VNPDAobj <- VNPDA(X.train,vy.train,X.test,rep(1,n.varbs),1,1,1,n.varbs,maxdepth=10)
  #Errors2[sim.rep] <- sum(abs(vy.test!=VNPDAobj$ClassPred))
  #RecordChoose2[sim.rep,] <- c(round(VNPDAobj$omega))
  #MCC2[sim.rep] <- Calc.MCC(predvec =RecordChoose2[sim.rep,], truevec = gamma.true)
  
  #Relabel data for to meet penLDA-FL input standards
  vy.tibs <- vy.train
  vy.tibs[vy.tibs==0] <- 2
  vy.new.tibs <- vy.test
  vy.new.tibs[vy.new.tibs==0] <- 2
  #Use training data to find optimal tunign parameters for penLDA-FL
  ValidError = matrix(0,nrow=length(Lambda_vec),ncol=length(Lambda2_vec))
  for(l in 1:length(Lambda_vec))
  {
    for(l2 in 1:length(Lambda2_vec)){
      
      out <- PenalizedLDA(X.train,vy.tibs,xte=X.train, lambda=Lambda_vec[l], K=1, type="ordered",lambda2 = Lambda2_vec[l2])
      ValidError[l,l2] <- sum(out$ypred!=vy.tibs)
      
    }
    
  }
  BestInd <- which(ValidError == min(ValidError), arr.ind = TRUE)
  BestLambda = Lambda_vec[BestInd[1,1]]
  BestLambda2 = Lambda2_vec[BestInd[1,2]]
  #Fit penLDA-FL with optimal tuning parameters
  out <- PenalizedLDA(X.train,vy.tibs,xte=X.test, lambda=BestLambda, K=1,type="ordered",lambda2 = BestLambda2)
  Errors3[sim.rep] <- sum(out$ypred != vy.new.tibs)
  RecordChoose3[sim.rep,which(abs(out$discrim) > quantile(abs(out$discrim),0.8))] <- 1
  MCC3[sim.rep] <- Calc.MCC(predvec =RecordChoose3[sim.rep,], truevec = gamma.true)
  
  #Fit randomForest
  rfObj <- randomForest(x=X.train,y=as.factor(vy.train))
  rf.pred <-predict(rfObj,newdata=X.test)
  rf2imp = importance(rfObj)
  Errors4[sim.rep] <- sum((as.numeric(rf.pred)-1)!=vy.test)
  RecordChoose4[sim.rep,which(c(rf2imp)!=0)] <- 1
  MCC4[sim.rep] <- Calc.MCC(predvec =RecordChoose4[sim.rep,], truevec = gamma.true)
  
  #Fit VLDA (run if you installed VaDA)
  #VLDAobj <- VLDA(vy=vy.train, mX=X.train, mXtest=X.test,r=0.98,kappa=(1/7))
  #Errors5[sim.rep] <- sum(abs(c(VLDAobj$y.hat)!=vy.test))
  #RecordChoose5[sim.rep,] <- c(VLDAobj$omega)
  #MCC5[sim.rep] <- Calc.MCC(predvec =RecordChoose5[sim.rep,], truevec = gamma.true)
  
  #Fit SparseLDA
  ValidErrorSDA = rep(0,6)
  for(l in 1:6)
  {
    outSDA <- sda(x=X.train,y=factor(vy.train), lambda=LambdaSparseLDA_vec[l], stop=2)
    predObj <- predict(outSDA, X.train)
    ValidErrorSDA[l] <- sum( (as.numeric(predObj$class) - 1) !=vy.train)
  }
  BestLambdaSDA = LambdaSparseLDA_vec[max(which(ValidErrorSDA==min(ValidErrorSDA)))]
  SDAobj <- sda(x=X.train,y=factor(vy.train), lambda=BestLambdaSDA, stop=2)
  predObj.SDA <- predict(SDAobj, X.test)
  Errors6[sim.rep] <- sum( (as.numeric(predObj.SDA$class) - 1) !=vy.test)
  RecordChoose6[sim.rep,SDAobj$varIndex] <- 1
  MCC6[sim.rep] <- Calc.MCC(predvec =RecordChoose6[sim.rep,], truevec = gamma.true)
  
  #Fit SVM-L2 with heurstically-determined cost parameter
  co <- heuristicC(X.train)
  SVMobj <- LiblineaR(data=X.train, target=vy.train, type=3, cost=co, bias=TRUE, verbose=FALSE)
  svmPred <- predict(SVMobj,X.test)
  Errors7[sim.rep] <- sum(svmPred$predictions!=vy.test)

  #Fit SVM-L1 with CV-optimised cost parameter
  SVMobj2 <- LiblineaR(data=X.train, target=vy.train, type=5, cost=0.4, bias=TRUE, verbose=FALSE)
  svmPred2 <- predict(SVMobj2,X.test)
  Errors8[sim.rep] <- sum(svmPred2$predictions!=vy.test)
  RecordChoose8[sim.rep,which(c(SVMobj2$W)!=0)] <- 1
  MCC8[sim.rep] <- Calc.MCC(predvec =RecordChoose8[sim.rep,], truevec = gamma.true)
  
  #Fit CNN
  ValidErrorCNN = rep(0,16)
  for(l in 1:16){
    
    CNN.fit.temp <- CNN.fit(Raw.X.train = X.train, Raw.Y.train = vy.train, Raw.X.test = X.train, filters.num = CNNsettings$FILTERS[l], kernels.num = CNNsettings$KERNELS[l])
    ValidErrorCNN <- sum(CNN.fit.temp$yhat!=vy.train)
    rm(CNN.fit.temp)
    
  }
  BestCNNset <- (which.min(ValidErrorCNN))[1]
  CNNobj <- CNN.fit(Raw.X.train = X.train, Raw.Y.train = vy.train, Raw.X.test = X.test, filters.num = CNNsettings$FILTERS[BestCNNset], kernels.num = CNNsettings$KERNELS[BestCNNset])
  Errors9[sim.rep] <- sum(CNNobj$yhat!=vy.test)
  
  #Fit GPDA-stat
  GPobjStat <- GPDA.sparse.Stat.Revised(mXtrain = X.train, vytrain = vy.train, mXtest = X.test, train.cycles = 5, test.cycles = 3, delta = 1)
  Errors10[sim.rep] <-sum(abs(round(c(GPobjStat$xi))!=vy.test))
  RecordChoose9[sim.rep,] <- round(c(GPobjStat$vw))
  MCC9[sim.rep] <- Calc.MCC(predvec =RecordChoose9[sim.rep,], truevec = gamma.true)
  
  cat("Completed repetition: ", sim.rep, "\n")
  
}



