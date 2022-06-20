BiocManager::install("ProData")
BiocManager::install("PROcess")
list.of.packages <- c("Matrix", "rpx", "MALDIquantForeign", "Rcpp", "RcppArmadillo","mvtnorm", "optimization", "psych", "nleqslv", "matrixStats", "pracma", "penalizedLDA", "e1071", "caret", "randomForest", "LiblineaR", "sparseLDA", "caret", "keras", "Rfast", "nloptr", "inline", "lbfgs", "cvTools", "mlbench", "sparseLDA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

library(ProData)
library(rpx)
library(MALDIquantForeign)
library(ggplot2)
library(randomForest)
library(plyr)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(optimization)
library(psych)
library(nleqslv)
library(matrixStats)
library(penalizedLDA)
library(pracma)
library(cvTools)
library(e1071)
library(randomForest)
library(mlbench)
library(PROcess)
library(caret)
library(keras)
library(matrixStats)
library(LiblineaR)
#library(VaDA)
library(sparseLDA)
library(MASS)

sourceCpp("GPDAnonStatFinal.cpp")
source("GPDARversionFinal.R")
source("CNNcode.R")

# Data is collected from labs in three countries
N = 114+167+20+10+17+34
y = c(rep(0,114),rep(1,167),rep(0,20), rep(1,10),rep(0,17),rep(1,34))
labind = c(rep(1,114+167),rep(2,20+10),rep(3,17+34))
table(y,labind)

# import data
i=1
mzf <- paste("N",paste(i),".mzML",sep="")
mzf <- file.path("Lab1/Neg", mzf)
ms <- importMzMl(mzf)
plot(ms[[1]])
for (i in 2:114){
  mzf <- paste("N",paste(i),".mzML",sep="")
  mzf <- file.path("Lab1/Neg", mzf)
  ms2 <- importMzMl(mzf)
  ms[[i]] = ms2[[1]]
}
for (i in 1:167){
  mzf <- paste("P",paste(i),".mzML",sep="")
  mzf <- file.path("Lab1/Pos", mzf)
  ms2 <- importMzMl(mzf)
  ms[[i+114]] = ms2[[1]]
}
for (i in 1:20){
  mzf <- paste("N",paste(i),".mzML",sep="")
  mzf <- file.path("Lab2/Neg", mzf)
  ms2 <- importMzMl(mzf)
  ms[[i+114+167]] = ms2[[1]]
}
for (i in 1:10){
  mzf <- paste("P",paste(i),".mzML",sep="")
  mzf <- file.path("Lab2/Pos", mzf)
  ms2 <- importMzMl(mzf)
  ms[[i+114+167+20]] = ms2[[1]]
}
for (i in 1:17){
  mzf <- paste("N",paste(i),".mzML",sep="")
  mzf <- file.path("Lab3/Neg", mzf)
  ms2 <- importMzMl(mzf)
  ms[[i+114+167+20+10]] = ms2[[1]]
}
for (i in 1:34){
  mzf <- paste("P",paste(i),".mzML",sep="")
  mzf <- file.path("Lab3/Pos", mzf)
  ms2 <- importMzMl(mzf)
  ms[[i+114+167+20+10+17]] = ms2[[1]]
}

# To process the data, follow the steps described in the paper and in this tutorial:
# https://cran.r-project.org/web/packages/MALDIquant/vignettes/MALDIquant-intro.pdf

# quality control
any(sapply(ms, isEmpty))
table(sapply(ms, length))

#trim to 3 to 15.5 kDA
tms = list()
for (i in 1:N){
  ind = (mass(ms[[i]])>=3000)&(mass(ms[[i]])<=15500)
  tms2 = createMassSpectrum(mass = mass(ms[[i]])[ind], intensity = intensity(ms[[i]])[ind], metaData = metaData(ms[[i]]))
  tms[[i]] = tms2
}
table(sapply(tms, length))

myfunc = function(object){
  sum(object@intensity <= 0L) 
}
table(sapply(tms, myfunc))

# variance stabalization
spectra <- transformIntensity(tms,method="sqrt")
table(sapply(spectra, myfunc))
## smoothing - NOTE: IT IS THE SMOOTHING FUNCTION WHICH CAUSES ZERO INTENSITIES!
spectra <- smoothIntensity(spectra, method="SavitzkyGolay",halfWindowSize=20)
table(sapply(spectra, myfunc))
# baseline subtraction
baseline <- estimateBaseline(spectra[[1]], method="SNIP",iterations=100)
plot(spectra[[1]])
lines(baseline, col="red", lwd=2)
spectra <- removeBaseline(tms, method="SNIP", iterations=100)
plot(spectra[[1]])
table(sapply(spectra, myfunc))
# calibration/normalization
spectra <- calibrateIntensity(spectra, method="TIC")
table(sapply(spectra, myfunc))

# peak detection
t = seq(3000,15500,.5)
peaks <- detectPeaks(spectra, method="MAD", halfWindowSize=10, SNR=2)
## warping/ALIGNMENT - NOTE: I DECREASED MINFREQUENCY TO 0.4 (so that some peaks may only be present in one group)
warpingFunctions <- determineWarpingFunctions(peaks, tolerance=0.003,minFrequency=0.4)
warpedSpectra <- warpMassSpectra(spectra, warpingFunctions)
warpedPeaks <- warpMassPeaks(peaks, warpingFunctions)
table(sapply(warpedSpectra, myfunc))
table(sapply(warpedSpectra, length))

### Processed Features in paper
fp = read.csv("41587_2020_644_MOESM3_ESM.csv")
TP = length(t)
detpeaks = colnames(fp)[2:(dim(fp)[2]-1)]
detpeak = as.numeric(substring(detpeaks,2))
peaksind = rep(0,TP)
for (i in 1:length(detpeak)){
  peaksind[sum(t<detpeak[i])] = 1
  peaksind[sum(t<detpeak[i])+1] = 1
}

source("intensityMatrix_v2.R")
featureMatrix <- intensityMatrix_v2(warpedSpectra,t)

## Plot data
TP = length(t)
N = length(y)
X = featureMatrix
X = sqrt(X) - mean(sqrt(X))


vy = as.numeric(as.factor(y[y<=1]))-1
X1 = X[vy==1,]
X0 = X[vy==0,]

combineXscale = 10*X


TRIALS <- 50
V<-5


n <- nrow(combineXscale)
p<-ncol(combineXscale)

Lambda_vec<-10^(seq(-3,2,1))
Lambda2_vec<-10^(seq(-3,2,1))
LambdaSparseLDA_vec<-c(1e-6,1e-4,1e-3,.1,1,10)
CNNsettings <- expand.grid(c(5,10,15,20), c(10,20,30,40))
names(CNNsettings) <- c("FILTERS","KERNELS")

Errors1 <- matrix(0,TRIALS,V)
Errors2 <- matrix(0,TRIALS,V)
Errors3 <- matrix(0,TRIALS,V)
Errors4 <- matrix(0,TRIALS,V)
Errors5 <- matrix(0,TRIALS,V)
Errors6 <- matrix(0,TRIALS,V)
Errors7 <- matrix(0,TRIALS,V)
Errors8 <- matrix(0,TRIALS,V)
Errors9 <- matrix(0,TRIALS,V)
Errors10 <- matrix(0,TRIALS,V)
Errors11 <- matrix(0,TRIALS,V)
Errors12 <- matrix(0,TRIALS,V)

RecordChoose1 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose2 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose3 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose4 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose5 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose6 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose8 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose9 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])


TPTN1 <- matrix(0, (TRIALS*V), 4)
TPTN2 <- matrix(0, (TRIALS*V), 4)
TPTN3 <- matrix(0, (TRIALS*V), 4)
TPTN4 <- matrix(0, (TRIALS*V), 4)
TPTN5 <- matrix(0, (TRIALS*V), 4)
TPTN6 <- matrix(0, (TRIALS*V), 4)
TPTN7 <- matrix(0, (TRIALS*V), 4)
TPTN8 <- matrix(0, (TRIALS*V), 4)
TPTN9 <- matrix(0, (TRIALS*V), 4)
TPTN10 <- matrix(0, (TRIALS*V), 4)
TPTN11 <- matrix(0, (TRIALS*V), 4)
TPTN12 <- matrix(0, (TRIALS*V), 4)

X1 = combineXscale[vy==1,]
X0 = combineXscale[vy==0,]

n <- length(vy)
ind_cv<-1
library(cvTools)


for (trial in 1:TRIALS) {
  
  set.seed(trial)
  a1 <- proc.time()[3]
  
  cvSets <- cvFolds(n, V)
  
  
  for(cv in 1:V){
    
    
    testInds <- cvSets$subsets[which(cvSets$which==cv)]
    trainInds <- (1:n)[-testInds]
    
    nTest <- length(testInds)
    nTrain <- length(trainInds)
    
    
    y.train <- vy[trainInds]
    X.train <- combineXscale[trainInds,]
    X.test <- combineXscale[testInds,]
    y.test <- vy[testInds]
    
    X.train1 <- X.train[y.train==1,]
    X.train0 <- X.train[y.train==0,]
    mean1 <- colMeans(X.train1)
    mean0 <- colMeans(X.train0)
    empirical.mZ <- t(apply(cbind(X.train,y.train),1,function(s) { ans=0; if(s[length(s)]==1){ ans = s[1:p] - mean1};  if(s[length(s)]==0){ ans = s[1:p] - mean0}; return(ans)  }))
    pvalues <- apply(X.train,2,function(s){ t.test(s[y.train==1],s[y.train==0],var.equal = TRUE)$p.value  })
    vw.init <- 1 - pvalues
    
    #Fit GPDA-Stat
    GPobjStat <- GPDA.sparse.Stat.Revised(mXtrain = X.train, vytrain = y.train, mXtest = X.test, train.cycles = 3, test.cycles = 3, delta = 1)
    Errors12[trial,cv] <- sum(abs(round(c(GPobjStat$xi))!=y.test))
    RecordChoose9[ind_cv,] <- round(c(GPobjStat$vw))
    TPTN12[ind_cv,] <- Calc.PNmatrix(round(c(GPobjStat$xi)), y.test)
    
    #Fit GPDA
    GPobj <- GPDA.sparse.NonStat.Proteomics(mXtrain = X.train, vytrain = y.train, mXtest = X.test, train.cycles = 5, test.cycles = 3, delta = 1)
    Errors1[trial,cv] <- sum(abs(round(c(GPobj$xi))!=y.test))
    RecordChoose1[ind_cv,] <- round(GPobj$vw)
    TPTN1[ind_cv,] <- Calc.PNmatrix(round(c(GPobj$xi)), y.test)

    #Fit VNPDA model(run if you installed VaDA)
    # VNPDAobj <- VNPDA(X.train,y.train,X.test,rep(1,p),1,1,1,p,maxdepth=10)
    # Errors2[trial,cv] <- sum(abs(y.test!=VNPDAobj$ClassPred))
    # RecordChoose2[ind_cv,] <- c(round(VNPDAobj$omega))
    # TPTN2[ind_cv,] <- Calc.PNmatrix(VNPDAobj$ClassPred, y.test)

    #Relabel data for to meet penLDA-FL input standards
    vy.tibs <- y.train
    vy.tibs[vy.tibs==0] <- 2
    vy.new.tibs <- y.test
    vy.new.tibs[vy.new.tibs==0] <- 2
    #Use training data to find optimal tuning parameters for penLDA-FL
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
    out <- PenalizedLDA(X.train,vy.tibs,xte=X.test, lambda=BestLambda, K=1,type="ordered",lambda2 = BestLambda2)
    Errors3[trial,cv] <- sum(out$ypred != vy.new.tibs)
    RecordChoose3[ind_cv,which(abs(out$discrim) > quantile(abs(out$discrim),0.8))] <- 1
    TPTN3[ind_cv,] <- Calc.PNmatrix(out$ypred, vy.new.tibs)

    #Fit randomForest
    rfObj <- randomForest(x=X.train,y=as.factor(y.train))
    rf.pred <-predict(rfObj,newdata=X.test)
    rf2imp = importance(rfObj)
    Errors4[trial,cv] <- sum((c(rf.pred)-1)!=y.test)
    RecordChoose4[ind_cv,which(c(rf2imp)!=0)] <- 1
    TPTN4[ind_cv,] <- Calc.PNmatrix((c(rf.pred)-1), y.test)

    #Fit VLDA (run if you installed VaDA)
    # VLDAobj <- VLDA(vy=y.train, mX=X.train, mXtest=X.test,r=0.98,kappa=(1/7))
    # Errors5[trial,cv] <- sum(abs(c(VLDAobj$y.hat)!=y.test))
    # RecordChoose5[ind_cv,] <- VLDAobj$omega
    # TPTN5[ind_cv,] <- Calc.PNmatrix(c(VLDAobj$y.hat), y.test)
    
    #Fit SparseLDA
    ValidErrorSDA = rep(0,6)
    for(l in 1:6)
    {
      outSDA <- sda(x=X.train,y=factor(y.train), lambda=LambdaSparseLDA_vec[l], stop=2)
      predObj <- predict(outSDA, X.train)
      ValidErrorSDA[l] <- sum( (as.numeric(predObj$class) - 1) !=y.train)
    }
    BestLambdaSDA = LambdaSparseLDA_vec[max(which(ValidErrorSDA==min(ValidErrorSDA)))]
    SDAobj <- sda(x=X.train,y=factor(y.train), lambda=BestLambdaSDA, stop=2)
    predObj.SDA <- predict(SDAobj, X.test)
    Errors6[trial,cv] <- sum( (as.numeric(predObj.SDA$class) - 1) !=y.test)
    TPTN6[ind_cv,] <- Calc.PNmatrix((as.numeric(predObj.SDA$class) - 1), y.test)
    RecordChoose6[ind_cv,SDAobj$varIndex] <- 1

    #Fit SVM-L2 with heurstically-determined cost parameter
    co <- heuristicC(X.train)
    SVMobj <- LiblineaR(data=X.train, target=y.train, type=3, cost=co, bias=TRUE, verbose=FALSE)
    svmPred <- predict(SVMobj,X.test)
    Errors7[trial,cv] <- sum(svmPred$predictions!=y.test)
    TPTN7[ind_cv,] <- Calc.PNmatrix(svmPred$predictions, y.test)
    RecordChoose6[ind_cv,which(abs(SVMobj$W[1:p])>quantile(abs(SVMobj$W[1:p]),0.7))] <- 1

    #Fit SVM-L1 with CV-optimised cost parameter
    co <- heuristicC(X.train)
    SVMobj2 <- LiblineaR(data=X.train, target=y.train, type=5, cost=0.4, bias=TRUE, verbose=FALSE)
    svmPred2 <- predict(SVMobj2,X.test)
    Errors8[trial,cv] <- sum(svmPred2$predictions!=y.test)
    TPTN8[ind_cv,] <- Calc.PNmatrix(svmPred2$predictions, y.test)
    RecordChoose8[ind_cv,which(c(SVMobj2$W[1:p])!=0)] <- 1

    #Fit Traditional LDA
    trainData <- data.frame(cbind(X.train[,peaksind==1],y.train))
    names(trainData) <- c(colnames(X.train)[peaksind==1],"yresponse")
    testData <- data.frame(cbind(X.test[,peaksind==1]))
    names(testData) <- colnames(X.train)[peaksind==1]
    testSVM.traditional <- testLDA.traditional <- lda(yresponse~.,data=trainData,)
    traditionalPred <- round(predict(testLDA.traditional,newdata=testData)$posterior[,2])
    TPTN9[ind_cv,] <- Calc.PNmatrix(traditionalPred, y.test)
    Errors9[trial,cv] <- sum(traditionalPred!=y.test)

    #Fit Traditional QDA
    QDApeaksind <- rep(0,p)
    QDApeaksind[sample(which(peaksind==1),100,replace=FALSE)] <-1
    QDAtrainData <- data.frame(cbind(X.train[,QDApeaksind==1],y.train))
    names(QDAtrainData) <- c(colnames(X.train)[QDApeaksind==1],"yresponse")
    QDAtestData <- data.frame(cbind(X.test[,QDApeaksind==1]))
    names(QDAtestData) <- colnames(X.train)[QDApeaksind==1]
    testQDA.traditional <- qda(yresponse~.,data=QDAtrainData,)
    traditionalQDAPred <- round(predict(testQDA.traditional,newdata=QDAtestData)$posterior[,2])
    TPTN10[ind_cv,] <- Calc.PNmatrix(traditionalQDAPred, y.test)
    Errors10[trial,cv] <- sum(traditionalQDAPred!=y.test)

    #Fit CNN
    ValidErrorCNN = rep(0,16)
    for(l in 1:16){

      CNN.fit.temp <- CNN.fit(Raw.X.train = X.train, Raw.Y.train = y.train, Raw.X.test = X.train, filters.num = CNNsettings$FILTERS[l], kernels.num = CNNsettings$KERNELS[l])
      ValidErrorCNN <- sum(CNN.fit.temp$yhat!=y.train)
      rm(CNN.fit.temp)

    }
    BestCNNset <- (which.min(ValidErrorCNN))[1]
    CNNobj <- CNN.fit(Raw.X.train = X.train, Raw.Y.train = y.train, Raw.X.test = X.test, filters.num = CNNsettings$FILTERS[BestCNNset], kernels.num = CNNsettings$KERNELS[BestCNNset])
    Errors11[trial,cv] <- sum(CNNobj$yhat!=y.test)
    TPTN11[ind_cv,] <- Calc.PNmatrix(CNNobj$yhat, y.test)
     rm(CNNobj)
     gc();gc()
      
    
    ind_cv <- ind_cv + 1
    #cat(trial,cv,Errors1[trial,cv],Errors2[trial,cv],Errors3[trial,cv],Errors4[trial,cv],Errors5[trial,cv],Errors6[trial,cv],Errors7[trial,cv],Errors8[trial,cv],nTest,"\n")
    
  }

  cat("Total errors for trial=", trial, "\n")
  cat("GPDA: ", sum(Errors1[trial,])/(N), "\n")
  cat("VNPDA: ", sum(Errors2[trial,])/(N), "\n")
  cat("penLDA-FL: ", sum(Errors3[trial,])/(N), "\n")
  cat("RandomForest: ", sum(Errors4[trial,])/(N), "\n")
  cat("VLDA: ", sum(Errors5[trial,])/(N), "\n")
  cat("SparseDA: ", sum(Errors6[trial,])/(N), "\n")
  cat("SVM-L2: ", sum(Errors7[trial,])/(N), "\n")
  cat("SVM-L1: ", sum(Errors8[trial,])/(N), "\n")
  cat("LDA-Traditional: ", sum(Errors9[trial,])/(N), "\n")
  cat("QDA-Traditional: ", sum(Errors10[trial,])/(N), "\n")
  cat("CNN: ", sum(Errors11[trial,])/(N), "\n")
  cat("GPDAstat: ", sum(Errors12[trial,])/(N), "\n")

}