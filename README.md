# Gaussian process discriminant analysis (GPDA)

Supplementary codes to: Non-stationary Gaussian process discriminant analysis with variable selection for high-dimensional functional data

* * *

# Prerequisite packages

Please ensure that you have the following R packages fully installed and functional before proceeding to run any code on this repository: Matrix; Rcpp; RcppArmadillo; mvtnorm; optimization; matrixStats; pracma; VaDA (link [here](https://github.com/weichangyu10/VaDA)); penalizedLDA; caret; randomForest; LiblineaR; sparseLDA; ggplot2

* * *

# Before running the codes... (Mac Users Only)
Please skip this step if you are not a Mac user. The package requires a Fortran compiler (through its RcppArmadillo dependency).
Chances are, your current Fortran compiler is not up-to-date. To update your Fortran compiler, simply follow the steps here: <br />
&nbsp;

1. In your Mac App Store, search "Xcode" and install. <br />
2. Open Terminal application. Type in

```{eval=FALSE}
xcode-select install
```
&nbsp; &nbsp;&nbsp;
and follow the instructions.<br />
&nbsp; &nbsp;&nbsp;
3. Click on the link [here](https://github.com/fxcoudert/gfortran-for-macOS/releases). Download the gfortan dmg file according to your MacOS version. <br />
&nbsp; &nbsp;&nbsp;
4. Open the dmg file, run the gfortran installer, follow all the instructions.

* * *

# Compute MCC and classification accuracy for simulation 1
```{r}
#Download "GPDARversionalFinal.R"
#Download "GPDAnonStatFinal.cpp"
#Download "GenSim.R"
#Download "RunSimulation.R"
#Change your simulation settings here. Possible values are 1, 2, 3, 4.
sim.set <- 1
source("RunSimulation.R")
```


* * *

# Plot results for simulation 1
```{r}
library(ggplot2)
#Simulation1
AccMATRIX1 <- 100 - cbind(Errors1,Errors2,Errors3,Errors4,Errors5,Errors6,Errors7,Errors8)/500*100
colnames(AccMATRIX1) <- c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L2pen","SVM-L1pen")
MCCMATRIX1 <- cbind(MCC1,MCC2,MCC3,MCC4,MCC5,MCC6,MCC8)*100
colnames(MCCMATRIX1) <- c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L1")
dat1 = reshape2::melt(AccMATRIX1,varnames=c("Iteration","Method"),value.name="ClassAccs")
names(dat1) <- c("Iteration","Classifier","Classification_Accs")
ClassError_plot1 <- ggplot(dat1, aes(x = Classifier, y = Classification_Accs,fill = Classifier)) + geom_boxplot()+ xlab(" ") + ylab("Classification Accuracy Percent") + theme_bw() + scale_fill_manual(values=c("yellow",rep("gray",7))) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "none", strip.text = element_text(size = 24))+theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5), strip.background = element_rect(color = "black", size = 1.5))+ggtitle("Simulation1")
ClassError_plot1

datMCC1 = reshape2::melt(MCCMATRIX1,varnames=c("Iteration","Method"),value.name="MCCval")
names(datMCC1) <- c("Iteration","Method","MCC")
MCC_plot1 <- ggplot(datMCC1, aes(x = Method, y = MCC,fill = Method)) + geom_boxplot()+ xlab(" ") + ylab("Matthews Correlation Coefficient") + theme_bw() + scale_fill_manual(values=c("yellow",rep("gray",6))) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "none", strip.text = element_text(size = 24))+theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5), strip.background = element_rect(color = "black", size = 1.5))+ggtitle("Simulation1")
MCC_plot1
```

* * *
# Bibliography
1. Yu, W., Wade, S., Bondell, H.D., Azizi, L. 2021. Non-stationary Gaussian process discriminant analysis with variable selection for high-dimensional functional data. <br />