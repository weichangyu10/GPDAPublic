# Gaussian process discriminant analysis (GPDA)

Vignette for: Non-stationary Gaussian process discriminant analysis with variable selection for high-dimensional functional data

* * *

# Prerequisite packages for all codes in this repository

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

# Toy example - Breast Cancer MS dataset
In this toy example, we will fit the GPDA model to the breast cancer MS dataset in the R package ProData. This dataset has been analyzed in section 5 of the reference paper. We begin with some pre-processing steps, which are similar to the steps described in Li et al. (2005). Note that outputs for the following example are embeded in Vignette.pdf.

Load required packages
```{r,results='hide',message=FALSE,tidy=TRUE, tidy.opts=list(width.cutoff=60)}
Package.Names <- c("Matrix","Rcpp","RcppArmadillo","mvtnorm","optimization","matrixStats","matrixStats","ggplot2","formatR")
lapply(Package.Names, library, character.only = TRUE)
sourceCpp("GPDAnonStatFinal.cpp")
source("GPDARversionFinal.R")
```

Load and import data
```{r,message=FALSE,tidy=TRUE, tidy.opts=list(width.cutoff=60)}
#BiocManager::install("PROcess")
#BiocManager::install("ProData")
library(ProData)
library(PROcess)
f45c <- system.file("f45c", package="ProData")
fs <- dir(f45c,full.names=TRUE)
data(f45cbmk)
SpecGrp <- pData(f45cbmk)
table(SpecGrp[,1])
```

Processed markers
```{r}
explevels = exprs(f45cbmk)
detpeaks = rownames(explevels)
```

Match spectra to sample data
```{r}
gi <- regexpr("i+[0-9]+", fs)
specName <- substr(fs, gi, gi + attr(gi, "match.length") - 1)
mt <- match(SpecGrp[, 2], toupper(specName))
```

Prepare classification response
```{r}
N = dim(SpecGrp)[1]
y = rep(0,N) # HER2 positive
y[SpecGrp[,1]=="B"] = 1 #healthy
y[SpecGrp[,1]=="C"] = 2 #ER/PR positive
y[SpecGrp[,1]=="D"] = 3 #single individual
```

Use the PROcess package for basic processing of the spectrum data. Process spectrum: baseline subtraction and renormalization
```{r}
spec <- rmBaseline(f45c,method="loess", bw=0.1)
spec <- spec[, mt] #match ordering to phenotype data
colnames(spec) <- SpecGrp[, 2]
prospec <- renorm(spec,cutoff = 1000)
prospec = prospec + 0.5
prospec = log(prospec)
```

Time points
```{r}
t2 = as.numeric(rownames(prospec))
TP2 = length(t2)
```

Plot processed spectrum data for healthy women
```{r}
df=data.frame(x = t(matrix(prospec[,y==1], 1,TP2*sum(y==1))), t = matrix(t2,TP2*sum(y==1),1), ind = c(1:sum(y==1))%x%rep(1,TP2))
df2=data.frame(x =apply(t(prospec[,y==1]),2,mean) , t = t2)
ggplot(df2, mapping = aes(x = t, y = x)) +
  geom_line(data = df, mapping = aes(x = t, y = x,  group = ind), colour = "GRAY", alpha = 1/2, size = 1/2) +
  geom_line(col= "red") +
  theme_bw() +
  xlab('m/z') +
  ylab('log intensity')
```

Plot processed spectrum data for HER2 positive
```{r}
df=data.frame(x = t(matrix(prospec[,y==0], 1,TP2*sum(y==0))), t = matrix(t2,TP2*sum(y==0),1), ind = c(1:sum(y==0))%x%rep(1,TP2))
df2=data.frame(x =apply(t(prospec[,y==0]),2,mean) , t = t2)
ggplot(df2, mapping = aes(x = t, y = x)) +
  geom_line(data = df, mapping = aes(x = t, y = x,  group = ind), colour = "GRAY", alpha = 1/2, size = 1/2) +
  geom_line(col= "red") +
  theme_bw() +
  xlab('m/z') +
  ylab('log intensity')
```

Plot mean spectrum across groups
```{r}
df=data.frame(x =c(apply(t(prospec[,y==0]),2,mean),apply(t(prospec[,y==1]),2,mean))  , t = rep(t2,2), Class = c(rep("HER2", TP2),rep("control", TP2)))
ggplot(df, mapping = aes(x = t, y = x)) +
  geom_line(data = df, mapping = aes(x = t, y = x, colour = Class))  +
  theme_bw() +
  xlab('m/z') +
  ylab('log intensity')
```

Subset to HER2 positive and healthy women only, then scale data
```{r}
vy = as.numeric(as.factor(y[y<=1]))-1
combineX = t(prospec[,y<=1])
combineXscale <- 10*(combineX - mean(combineX))
p=ncol(combineXscale)
n <- length(vy)
```
    
Fit GPDA model. Description of the input arguments: mXtrain is a n.train by p matrix of predictors for training dataset. vytrain is a vector of size n.train of class labels for training dataset. mXtest is an n.test by p matrix of predictors for testing dataset. train.cycles is the number of VB cycles for posterior inference phase of the algorithm. test.cycles is the number of VB cycles for classification phase of the algorithm. delta is the distance between adjacent locations.
```{r}
GPobj <- GPDA.sparse.NonStat(mXtrain = combineXscale, vytrain = vy, mXtest = combineXscale, train.cycles = 5, test.cycles = 3, delta = 1)
```

Obtain predicted class labels
```{r}
#c(GPobj$xi)
```

Obtain posterior probability of gamma = 1 for all locations
```{r}
#c(GPobj$vw)
```

Plot posterior expectation of location-varying length scale of first observation.
```{r}
#Observation-specific log length scale perturbation
GPobj$zeta
#Common component of log length scale perturbation
#c(GPobj$mS)
#Posterior expectation of length-scale for first observation
LengthScaleObs1 <- exp(c(GPobj$mS)+GPobj$zeta[1])

df=data.frame(x = LengthScaleObs1, t = t2[1:(TP2-1)])
ggplot(df, mapping = aes(x = t, y = x)) +
  geom_line(data = df, mapping = aes(x = t, y = x))  +
  theme_bw() +
  xlab('m/z') +
  ylab('Length scale')
```

Plot posterior expectation of the first observation-specific latent process.
```{r}
#Diagonal entries of \Sigma_{z_1}
#GPobj$Sigma_Z[,1,1]
#m_{z_1}
#GPobj$mZ[1,]
z <- GPobj$mZ[1,]+0.5*GPobj$Sigma_Z[,1,1]
df=data.frame(x = z, t = t2)
ggplot(df, mapping = aes(x = t, y = x)) +
  geom_line(data = df, mapping = aes(x = t, y = x))  +
  theme_bw() +
  xlab('m/z') +
  ylab('E_q[Z(t)]')
```

* * *
# Bibliography
1. Yu, W., Wade, S., Bondell, H.D., Azizi, L. 2021. Non-stationary Gaussian process discriminant analysis with variable selection for high-dimensional functional data. <br />
