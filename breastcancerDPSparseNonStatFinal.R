library(ProData)
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
library(LiblineaR)
library(VaDA)
library(sparseLDA)
library(MASS)
sourceCpp("GPDAnonStatFinal.cpp")
source("GPDARversionFinal.R")

# Import raw spectrum data
f45c <- system.file("f45c", package="ProData")
fs <- dir(f45c,full.names=TRUE)

# Import sample data
data(f45cbmk)
SpecGrp <- pData(f45cbmk)
table(SpecGrp[,1])
#55 are HER2 positive (A), 64 are normal healthy women (B), 35 are ER/PR positive (mostly) (C) 
#and 13 samples are from a single healthy woman (D)

# Processed markers
explevels = exprs(f45cbmk)
detpeaks = rownames(explevels)

# Match spectra to sample data
gi <- regexpr("i+[0-9]+", fs)
specName <- substr(fs, gi, gi + attr(gi, "match.length") - 1)
mt <- match(SpecGrp[, 2], toupper(specName))

# classification response
N = dim(SpecGrp)[1]
y = rep(0,N) # HER2 positive
y[SpecGrp[,1]=="B"] = 1 #healthy
y[SpecGrp[,1]=="C"] = 2 #ER/PR positive
y[SpecGrp[,1]=="D"] = 3 #single individual

##################################################
## Use the PROcess package for basic processing of the specturm data
# process spectrum: baseline subtraction and renormalized
spec <- rmBaseline(f45c,method="loess", bw=0.1)
spec <- spec[, mt] #match ordering to phenotype data
colnames(spec) <- SpecGrp[, 2]
prospec <- renorm(spec,cutoff = 1000)
## log transformation
min(min(prospec))
prospec = prospec + 0.5
min(min(prospec))
prospec = log(prospec)

# time points
t2 = as.numeric(rownames(prospec))
TP2 = length(t2)

detpeak = as.numeric(substring(detpeaks,2))
peaksind = rep(0,TP2)
for (i in 1:length(detpeak)){
  peaksind[sum(t2<detpeak[i])] = 1
  peaksind[sum(t2<detpeak[i])+1] = 1
}

vy = as.numeric(as.factor(y[y<=1]))-1
combineX = t(prospec[,y<=1])
p=ncol(combineX);


n <- length(vy)

TRIALS = 50
V = 5

combineXscale <- 10*(combineX - mean(combineX))

Lambda_vec<-10^(seq(-3,2,1))
Lambda2_vec<-10^(seq(-3,2,1))
LambdaSparseLDA_vec<-c(1e-6,1e-4,1e-3,.1,1,10)


X1 = combineXscale[vy==1,]
X0 = combineXscale[vy==0,]
normality.p1 <- apply(X1,2,function(s){ shapiro.test(s)$p.value })
normality.p0 <- apply(X0,2,function(s){ shapiro.test(s)$p.value })
mean(normality.p1<0.05)
mean(normality.p0<0.05)



p <- ncol(X1)
ind_cv<-1
library(cvTools)

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

RecordChoose1 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose2 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose3 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose4 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose5 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose6 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])
RecordChoose8 <- matrix(0, (TRIALS*V), dim(combineXscale)[2])

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

X1 = combineXscale[vy==1,]
X0 = combineXscale[vy==0,]

n <- length(vy)
ind_cv<-1

for(trial in 1:TRIALS) {
  
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
    
    GPobj <- GPDA.sparse.NonStat(mXtrain = X.train, vytrain = y.train, mXtest = X.test, train.cycles = 5, test.cycles = 3, delta = 1)
    Errors1[trial,cv] <- sum(abs(round(c(GPobj$xi))!=y.test))
    RecordChoose1[ind_cv,] <- round(GPobj$vw)
    TPTN1[ind_cv,] <- Calc.PNmatrix(round(c(GPobj$xi)), y.test)


    VNPDAobj <- VNPDA(X.train,y.train,X.test,rep(1,p),1,1,1,p,maxdepth=10)
    Errors2[trial,cv] <- sum(abs(y.test!=VNPDAobj$ClassPred))
    RecordChoose2[ind_cv,] <- c(round(VNPDAobj$omega))
    TPTN2[ind_cv,] <- Calc.PNmatrix(VNPDAobj$ClassPred, y.test)

    vy.tibs <- y.train
    vy.tibs[vy.tibs==0] <- 2
    vy.new.tibs <- y.test
    vy.new.tibs[vy.new.tibs==0] <- 2

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

    rfObj <- randomForest(x=X.train,y=as.factor(y.train))
    rf.pred <-predict(rfObj,newdata=X.test)
    rf2imp = importance(rfObj)
    Errors4[trial,cv] <- sum((c(rf.pred)-1)!=y.test)
    RecordChoose4[ind_cv,which(c(rf2imp)!=0)] <- 1
    TPTN4[ind_cv,] <- Calc.PNmatrix((c(rf.pred)-1), y.test)

    VLDAobj <- VLDA(vy=y.train, mX=X.train, mXtest=X.test,r=0.98,kappa=(1/7))
    Errors5[trial,cv] <- sum(abs(c(VLDAobj$y.hat)!=y.test))
    RecordChoose5[ind_cv,] <- VLDAobj$omega
    TPTN5[ind_cv,] <- Calc.PNmatrix(c(VLDAobj$y.hat), y.test)

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
    
    co <- heuristicC(X.train)
    SVMobj <- LiblineaR(data=X.train, target=y.train, type=3, cost=co, bias=TRUE, verbose=FALSE)
    svmPred <- predict(SVMobj,X.test)
    Errors7[trial,cv] <- sum(svmPred$predictions!=y.test)
    TPTN7[ind_cv,] <- Calc.PNmatrix(svmPred$predictions, y.test)

    SVMobj2 <- LiblineaR(data=X.train, target=y.train, type=5, cost=0.4, bias=TRUE, verbose=FALSE)
    svmPred2 <- predict(SVMobj2,X.test)
    Errors8[trial,cv] <- sum(svmPred2$predictions!=y.test)
    TPTN8[ind_cv,] <- Calc.PNmatrix(svmPred2$predictions, y.test)
    RecordChoose8[ind_cv,which(c(SVMobj2$W[1:p])!=0)] <- 1

    trainData <- data.frame(cbind(X.train[,peaksind==1],y.train))
    names(trainData) <- c(colnames(X.train)[peaksind==1],"yresponse")
    testData <- data.frame(cbind(X.test[,peaksind==1]))
    names(testData) <- colnames(X.train)[peaksind==1]
    testLDA.traditional <- lda(yresponse~.,data=trainData,)
    traditionalPred <- round(predict(testLDA.traditional,newdata=testData)$posterior[,2])
    TPTN9[ind_cv,] <- Calc.PNmatrix(traditionalPred, y.test)
    Errors9[trial,cv] <- sum(traditionalPred!=y.test)

    QDApeaksind <- rep(0,p)
     QDApeaksind[sample(which(peaksind==1),30,replace=FALSE)] <-1
     QDAtrainData <- data.frame(cbind(X.train[,QDApeaksind==1],y.train))
     names(QDAtrainData) <- c(colnames(X.train)[QDApeaksind==1],"yresponse")
     QDAtestData <- data.frame(cbind(X.test[,QDApeaksind==1]))
     names(QDAtestData) <- colnames(X.train)[QDApeaksind==1]
    testQDA.traditional <- qda(yresponse~.,data=QDAtrainData,)
    traditionalQDAPred <- round(predict(testQDA.traditional,newdata=QDAtestData)$posterior[,2])
    TPTN10[ind_cv,] <- Calc.PNmatrix(traditionalQDAPred, y.test)
    Errors10[trial,cv] <- sum(traditionalQDAPred!=y.test)
    
    ind_cv <- ind_cv + 1
    
  }
  cat("Total errors for trial=", trial, "\n")
  cat("GPDA: ", sum(Errors1[trial,])/(n), "\n")
  cat("VNPDA: ", sum(Errors2[trial,])/(n), "\n")
  cat("penLDA-FL: ", sum(Errors3[trial,])/(n), "\n")
  cat("RandomForest: ", sum(Errors4[trial,])/(n), "\n")
  cat("VLDA: ", sum(Errors5[trial,])/(n), "\n")
  cat("SparseDA: ", sum(Errors6[trial,])/(n), "\n")
  cat("SVM-L2: ", sum(Errors7[trial,])/(n), "\n")
  cat("SVM-L1: ", sum(Errors8[trial,])/(n), "\n")
  cat("LDA-Traditional: ", sum(Errors9[trial,])/(n), "\n")
  cat("QDA-Traditional: ", sum(Errors10[trial,])/(n), "\n")
  
  
}


TP1.transformed <- rep(0,TRIALS)
TP2.transformed <- rep(0,TRIALS)
TP3.transformed <- rep(0,TRIALS)
TP4.transformed <- rep(0,TRIALS)
TP5.transformed <- rep(0,TRIALS)
TP6.transformed <- rep(0,TRIALS)
TP7.transformed <- rep(0,TRIALS)
TP8.transformed <- rep(0,TRIALS)
TP9.transformed <- rep(0,TRIALS)
TP10.transformed <- rep(0,TRIALS)
TN1.transformed <- rep(0,TRIALS)
TN2.transformed <- rep(0,TRIALS)
TN3.transformed <- rep(0,TRIALS)
TN4.transformed <- rep(0,TRIALS)
TN5.transformed <- rep(0,TRIALS)
TN6.transformed <- rep(0,TRIALS)
TN7.transformed <- rep(0,TRIALS)
TN8.transformed <- rep(0,TRIALS)
TN9.transformed <- rep(0,TRIALS)
TN10.transformed <- rep(0,TRIALS)
for(trial in 1:TRIALS){
  
  startid <- (trial-1)*5+1
  
  tempTP1 <- colSums(TPTN1[startid:(startid+4),])
  TP1.transformed[trial] <-  tempTP1[1]/(tempTP1[1]+tempTP1[3])
  TN1.transformed[trial] <-  tempTP1[2]/(tempTP1[2]+tempTP1[4])
  
  tempTP2 <- colSums(TPTN2[startid:(startid+4),])
  TP2.transformed[trial] <-  tempTP2[1]/(tempTP2[1]+tempTP2[3])
  TN2.transformed[trial] <-  tempTP2[2]/(tempTP2[2]+tempTP2[4])
  
  tempTP3 <- colSums(TPTN3[startid:(startid+4),])
  TP3.transformed[trial] <-  tempTP3[1]/(tempTP3[1]+tempTP3[3])
  TN3.transformed[trial] <-  tempTP3[2]/(tempTP3[2]+tempTP3[4])
  
  tempTP4 <- colSums(TPTN4[startid:(startid+4),])
  TP4.transformed[trial] <-  tempTP4[1]/(tempTP4[1]+tempTP4[3])
  TN4.transformed[trial] <-  tempTP4[2]/(tempTP4[2]+tempTP4[4])
  
  tempTP5 <- colSums(TPTN5[startid:(startid+4),])
  TP5.transformed[trial] <-  tempTP5[1]/(tempTP5[1]+tempTP5[3])
  TN5.transformed[trial] <-  tempTP5[2]/(tempTP5[2]+tempTP5[4])
  
  tempTP6 <- colSums(TPTN6[startid:(startid+4),])
  TP6.transformed[trial] <-  tempTP6[1]/(tempTP6[1]+tempTP6[3])
  TN6.transformed[trial] <-  tempTP6[2]/(tempTP6[2]+tempTP6[4])
  
  tempTP7 <- colSums(TPTN7[startid:(startid+4),])
  TP7.transformed[trial] <-  tempTP7[1]/(tempTP7[1]+tempTP7[3])
  TN7.transformed[trial] <-  tempTP7[2]/(tempTP7[2]+tempTP7[4])
  
  tempTP8 <- colSums(TPTN8[startid:(startid+4),])
  TP8.transformed[trial] <-  tempTP8[1]/(tempTP8[1]+tempTP8[3])
  TN8.transformed[trial] <-  tempTP8[2]/(tempTP8[2]+tempTP8[4])
  
  tempTP9 <- colSums(TPTN9[startid:(startid+4),])
  TP9.transformed[trial] <-  tempTP9[1]/(tempTP9[1]+tempTP9[3])
  TN9.transformed[trial] <-  tempTP9[2]/(tempTP9[2]+tempTP9[4])
  
  tempTP10 <- colSums(TPTN10[startid:(startid+4),])
  TP10.transformed[trial] <-  tempTP10[1]/(tempTP10[1]+tempTP10[3])
  TN10.transformed[trial] <-  tempTP10[2]/(tempTP10[2]+tempTP10[4])
}

if(FALSE){
  
#Compute true positive and true negative rates
TPmatrix <- cbind(TP1.transformed,TP2.transformed,TP3.transformed,TP4.transformed,TP5.transformed,TP6.transformed,TP7.transformed,TP8.transformed,TP9.transformed,TP10.transformed)
TPmatrix <- c(TPmatrix)
TPmatrix <- matrix(TPmatrix,nrow=50,ncol=10)
colnames(TPmatrix) <- c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L2","SVM-L1","LDA-Traditional","QDA-Traditional")
TPmatrix <- 100*TPmatrix
datTP = reshape2::melt(TPmatrix,varnames=c("Iteration","Method"),value.name="FalsePositiveRate")
names(datTP) <- c("Iteration","Classifier","TruePositiveRate")
TP_plot <- ggplot(datTP, aes(x = Classifier, y = TruePositiveRate,fill = Classifier)) + geom_boxplot()+ xlab(" ") + ylab("True positive rate (%)") + theme_bw() + scale_fill_manual(values=c("yellow",rep("gray",9))) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "none", strip.text = element_text(size = 24))+theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5), strip.background = element_rect(color = "black", size = 1.5))+ggtitle("Breast cancer")
TP_plot

TNmatrix <- cbind(TN1.transformed,TN2.transformed,TN3.transformed,TN4.transformed,TN5.transformed,TN6.transformed,TN7.transformed,TN8.transformed,TN9.transformed,TN10.transformed)
TNmatrix <- c(TNmatrix)
TNmatrix <- matrix(TNmatrix,nrow=50,ncol=10)
colnames(TNmatrix) <- c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L2","SVM-L1","LDA-Traditional","QDA-Traditional")
TNmatrix <- 100*TNmatrix
datTN = reshape2::melt(TNmatrix,varnames=c("Iteration","Method"),value.name="TrueNegativeRate")
names(datTN) <- c("Iteration","Classifier","TrueNegativeRate")
TN_plot <- ggplot(datTN, aes(x = Classifier, y = TrueNegativeRate,fill = Classifier)) + geom_boxplot()+ xlab(" ") + ylab("True negative rate (%)") + theme_bw() + scale_fill_manual(values=c("yellow",rep("gray",9))) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "none", strip.text = element_text(size = 24))+theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5), strip.background = element_rect(color = "black", size = 1.5))+ggtitle("Breast cancer")
TN_plot


#Compute time point selection rates
GPDA.selectRate <- colMeans(RecordChoose1)*100
VNPDA.selectRate <- colMeans(RecordChoose2)*100
penLDA.selectRate <- colMeans(RecordChoose3)*100
rf.selectRate <- colMeans(RecordChoose4)*100
VLDA.selectRate <- colMeans(RecordChoose5)*100
SparseLDA.selectRate <- colMeans(RecordChoose6)*100
SVM.selectRate <- colMeans(RecordChoose8)*100
Traditional.selectRate <- peaksind*100

dat.varb.select <- cbind(GPDA.selectRate,VNPDA.selectRate,penLDA.selectRate,rf.selectRate,VLDA.selectRate,SparseLDA.selectRate,SVM.selectRate, Traditional.selectRate)
datVarbSelect.breast = reshape2::melt(dat.varb.select,varnames=c("Iteration","Method"),value.name="SelectionRate")
datVarbSelect.breast[,2] <- rep(c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L1","LDA-Traditional"),rep(p,8))
datVarbSelect.breast[,2] <- factor(datVarbSelect.breast[,2],levels=c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L1","LDA-Traditional"))
names(datVarbSelect.breast) <- c("Timepoint","Method","SelectionRate")
VarbSelect.HeatMap.breast <- ggplot(datVarbSelect.breast, aes(Method, Timepoint, fill= SelectionRate)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "bottom", legend.text=element_text(size=12)) +ggtitle("Breast cancer") + ylab("m/Z location")+ scale_y_continuous(breaks=seq(1000,20000,2000))
VarbSelect.HeatMap


#Compute classification accuracy rates
Accmatrix <- 100 - cbind(rowSums(Errors1),rowSums(Errors2),rowSums(Errors3), rowSums(Errors4),rowSums(Errors5), rowSums(Errors6), rowSums(Errors7) , rowSums(Errors8) , rowSums(Errors9) , rowSums(Errors10))/n*100
colnames(Accmatrix) <- c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L2","SVM-L1","LDA-Traditional","QDA-Traditional")
datAcc = reshape2::melt(Accmatrix,varnames=c("Iteration","Method"),value.name="ClassificationAccuracy")
names(datAcc) <- c("Iteration","Classifier","ClassificationAccuracy")
Acc_plot <- ggplot(datAcc, aes(x = Classifier, y = ClassificationAccuracy,fill = Classifier)) + geom_boxplot()+ xlab(" ") + ylab("Classification accuracy (%)") + theme_bw() + scale_fill_manual(values=c("yellow",rep("gray",9))) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "none", strip.text = element_text(size = 24))+theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5), strip.background = element_rect(color = "black", size = 1.5))+ggtitle("Breast cancer")
Acc_plot


#Expected log length scale of first observation
LogLengthScale = GPobj$mS + GPobj$zeta[1]
dataframeLogLengthScale <- data.frame(LogLengthScale,1:(p-1))
names(dataframeLogLengthScale) <- c("value","Timepoint")
ExpectedNuplots <- ggplot(dataframeLogLengthScale, aes(x = Timepoint, y = value)) + geom_line(size = 0.5) +scale_color_manual(values = c("red")) + theme_bw() + ggtitle("Expected log length scale of first observation") + theme(legend.position = "bottom", strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))

#Expected Z of first observation
EZ1 = GPobj$mZ[1,]
dataframeZ <- data.frame(EZ1,1:p)
names(dataframeZ) <- c("value","Timepoint")
ExpectedZplots <- ggplot(dataframeZ, aes(x = Timepoint, y = value)) + geom_line(size = 0.5) +scale_color_manual(values = c("red")) + theme_bw() + ggtitle("E_q(Z_1)") + theme(legend.position = "bottom", strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))


#Do the methods pick up the genes selected by traditional peak detection?
c(mean(colMeans(RecordChoose1)[peaksind==1]),mean(colMeans(RecordChoose2)[peaksind==1]), mean(colMeans(RecordChoose3)[peaksind==1]), mean(colMeans(RecordChoose4)[peaksind==1]), mean(colMeans(RecordChoose5)[peaksind==1]), mean(colMeans(RecordChoose6)[peaksind==1]), mean(colMeans(RecordChoose8)[peaksind==1]) )
c(sd(colMeans(RecordChoose1)[peaksind==1]),sd(colMeans(RecordChoose2)[peaksind==1]), sd(colMeans(RecordChoose3)[peaksind==1]), sd(colMeans(RecordChoose4)[peaksind==1]), sd(colMeans(RecordChoose5)[peaksind==1]), sd(colMeans(RecordChoose6)[peaksind==1]), sd(colMeans(RecordChoose8)[peaksind==1]) )

}