library(Rfast)
library(nloptr)
library(inline)
library(lbfgs)

FitEpsilonRversion <- function(X,y,mZ,m_mu1,m_mu0,m_mu,VarMu1,VarMu0,VarMu,VarZ,vw,Aeps = 0.00001,Beps = 0.000001){

  n <- nrow(mZ)
  n1 <- sum(y)
  n0 <- n - n1
  p <- ncol(mZ)
  
  ResidualDataSq1.vec <- (c(X[y==1,]) - c(mZ[y==1,]) - rep(m_mu1,rep(n1,p)))^2 + rep(VarMu1,rep(n1,p)) + c(VarZ[y==1,])
  ResidualDataSq0.vec <- (c(X[y==0,]) - c(mZ[y==0,]) - rep(m_mu0,rep(n0,p)))^2 + rep(VarMu0,rep(n0,p)) + c(VarZ[y==0,])
  ResidualDataSq.vec <- (c(X) - c(mZ) - rep(m_mu,rep(n,p)))^2 + rep(VarMu,rep(n,p)) + c(VarZ)
  
  Seps1 <- Beps + 0.5*vw*group(ResidualDataSq1.vec, rep(1:p,rep(n1,p)),"sum")
  Seps0 <- Beps + 0.5*vw*group(ResidualDataSq0.vec, rep(1:p,rep(n0,p)),"sum")
  Seps <- Beps + 0.5*(1-vw)*group(ResidualDataSq.vec, rep(1:p,rep(n,p)),"sum")
  
  Reps1 <- Aeps + 0.5*n1*vw
  Reps0 <- Aeps + 0.5*n0*vw
  Reps <- Aeps + 0.5*n*(1-vw)
  
  EDepsInv <- cbind(Reps1/Seps1,Reps0/Seps0,Reps/Seps)
  ELogeps <- cbind(log(Seps1) - digamma(Reps1),log(Seps0) - digamma(Reps0),log(Seps) - digamma(Reps))
  
  return(list(EDepsInv=EDepsInv, ELogeps = ELogeps))
  
}

OptimizeLambda <- function(ENusqTotal, ENuNuTotal, m_tau2_inv, mu.lambda, s2.lambda, samp.size, delta){
  
  p <- length(ENusqTotal) + 1
  
  coef1 <- sum(ENusqTotal[1:(p-2)])
  coef2 <- sum(ENuNuTotal)
  coef3 <- sum(ENusqTotal[2:(p-1)])
  
  fr <- function(x){
    
    lambdatilde <- x
    c0sq <- exp(lambdatilde)/(2*delta)
    c1sq <- c0sq + 0.5*delta*exp(-lambdatilde) - 1
    c0c1 <- 0.5*(1- exp(lambdatilde)/delta)
    OBJ <- 0.5*samp.size*(p-2)*lambdatilde - 0.5*m_tau2_inv*(c1sq*coef1 + 2*c0c1*coef2 + c0sq*coef3) - (lambdatilde-mu.lambda)^2/(2*s2.lambda)
    return(-OBJ)
  }
  
  grr <- function(x){
    
    lambdatilde <- x
    dc0sq <- exp(lambdatilde)/(2*delta)
    dc1sq <- dc0sq - 0.5*delta*exp(-lambdatilde)
    dc0c1 <- -dc0sq
    GRAD.OBJ <- 0.5*samp.size*(p-2) - 0.5*m_tau2_inv*(dc1sq*coef1 + 2*dc0c1*coef2 + dc0sq*coef3) - (lambdatilde-mu.lambda)/s2.lambda
    return(-GRAD.OBJ)
    
  }
  x.guess <- rnorm(1,sd=5)
  optimObj <- optim(x.guess, fr, grr, method = "L-BFGS-B",
                    lower = rep(-500, p), upper = rep(500, p), control = list(maxit = 200000))
  return(optimObj$par)
}

OptBeta <- function(sumW, sumAdjW, mu.beta, sigma2.beta, mu.alpha, sigma2.alpha, p){
  
  r.init <- c(-log(p), 0.3)
  OptimFunc <- function(r){
    
    alpha <- r[1]
    beta <- r[2]
    K <- beta/4
    L <- 0.5*(beta - alpha)
    eta1 <- exp(K)*cosh(L) + sqrt(exp(2*K)*(cosh(L))^2 - 2*sinh(2*K))
    eta2 <- exp(K)*cosh(L) - sqrt(exp(2*K)*(cosh(L))^2 - 2*sinh(2*K))
    logZ <- p*log(eta1) + log(1 + (eta2/eta1)^p)
    
    ObjecValue <- beta*sumAdjW - alpha*sumW + (K-L)*p - logZ - (log(beta) - mu.beta)^2/(2*sigma2.beta) - log(beta) - (alpha - mu.alpha)^2/(2*sigma2.alpha)
    
    return(ObjecValue)
  }
  
  FinalAns <- optim_sa(OptimFunc,start = r.init,maximization = TRUE,lower = c(-1000,0),upper=c(1000,50), control = list(nlimit=1000))
  
  return(FinalAns$par)
}

OptimizeZeta3 <- function(mZ, SigmaZ, mS, VarS, sigma2.zeta, m_tau_inv, delta){
  
  p <- ncol(mZ)
  n <- nrow(mZ)
  
  Ezsq <- (mZ^2 + t(SigmaZ[,1,]))
  Ezz <- (mZ[,1:(p-1)]*mZ[,2:p] + t(SigmaZ[1:(p-1),2,]))
  
  ES <- mS+0.5*VarS
  ENegS <- -mS+0.5*VarS
  
  x.guess <- rt(n,df=5)
  
  ans <- rep(0,n)
  
  minAns <- log(delta) - min(mS)
  
  for(i in 1:n){
    
    fr <- function(x){
      
      zeta.i <- x
      repzeta <- rep(zeta.i,p-1)
      b0sqEZsq1 <- 1/(2*delta)*sum( exp( repzeta + ES + log(Ezsq[i,2:p])) )
      b0sqEZsq2 <- 1/(2*delta)*sum( exp( repzeta + ES + log(Ezsq[i,1:(p-1)])) )
      b1sqEZsq <- b0sqEZsq2 + 0.5*delta*sum(exp(-repzeta + ENegS + log(Ezsq[i,1:(p-1)]))) - sum(Ezsq[i,1:(p-1)])
      b0b1Ezz <- 0.5*sum(Ezz[i,]) - 1/(2*delta)*sum(exp(repzeta + ES)*Ezz[i,])
      
      OBJ <- 0.5*(p-1)*zeta.i - 0.5*m_tau_inv*( b1sqEZsq + 2*b0b1Ezz + b0sqEZsq1 ) - 0.5/(sigma2.zeta)*zeta.i^2
      return(OBJ)
    }
    
    grr <- function(x){
      
      zeta.i <- x
      repzeta <- rep(zeta.i,p-1)
      db0sqEZsq1 <- 1/(2*delta)*sum(exp( repzeta + ES + log(Ezsq[i,2:p])))
      db0sqEZsq2 <- 1/(2*delta)*sum(exp( repzeta + ES + log(Ezsq[i,1:(p-1)])))
      db1sqEZsq <- db0sqEZsq2 - 0.5*delta*sum(exp( -repzeta  + ES + log(Ezsq[i,1:(p-1)]) ))
      db0b1Ezz <- - 1/(2*delta)*sum(exp( repzeta + ES)*Ezz[i,])
      gradCurr <- 0.5*(p-1) - 0.5*m_tau_inv*(db1sqEZsq + 2*db0b1Ezz + db0sqEZsq1) - zeta.i/sigma2.zeta
      return(gradCurr)
      
    }
    
    optimObj <- optim(x.guess[i], fr, grr, method = "BFGS", control = list(maxit = 50, fnscale = -1))
    
    ans[i] <- optimObj$par
    if(ans[i] < minAns + 0.01){
      
      ans[i] <- minAns + 0.01
      
    }
    
  }
  
  return(ans)
}

SetInitialXi <- function(mXtrain, vytrain, mXtest ,alpha=0.90, delta=0){

  initial.y.dist <- table(vytrain)/length(vytrain)
  
  trainMeans <- apply(mXtrain, 1, mean)
  mXtrain <- mXtrain-trainMeans
  Pred.y <- array(vytrain, dim=c(length(alpha), length(delta), length(vytrain)))
  mXtest <- mXtest-trainMeans
  
  Pred.y.new <- array(0, dim=c(length(alpha), length(delta),ncol(mXtest)))
  initialXival <- array(0, dim=c(length(alpha), length(delta), ncol(mXtest), length(initial.y.dist)))
  
  Y <- model.matrix(~factor(vytrain)-1)
  xbar <- scale(mXtrain%*%Y, FALSE, table(vytrain)) 
  mXtrainstar <- mXtrain-xbar[, unclass(factor(vytrain))]
  
  XTX <- t(mXtrainstar)%*%mXtrainstar  
  SVBobj <- svd(XTX)
  SVDval <- SVBobj$d                     
  SVDvalPos <- seq(SVDval)[SVDval > 0] 
  SVDval <- sqrt(SVDval[SVDvalPos])
  
  Vmat <- scale(mXtrainstar%*%SVBobj$u[, SVDvalPos], FALSE, SVDval)
  SVDval <- SVDval/sqrt(length(vytrain))
  Vmattxbar <- t(Vmat)%*%xbar
  FinalErr <- matrix(0, length(alpha), length(delta))
  
  for(i in seq(along=alpha)){
    tempObj <- Vmattxbar*(1/(SVDval^2*alpha[i]+(1-alpha[i]))-1/(1-alpha[i]))
    CM <- Vmat%*%tempObj+xbar/(1-alpha[i])
    for(j in seq(along=delta)){
      CM1 <- I(abs(CM)>delta[j])*(abs(CM)-delta[j])*sign(CM)
      
      yprime <- t(Vmat)%*%CM1
      intvec <- -(alpha[i]*apply(yprime^2*SVDval^2, 2, sum)+(1-alpha[i])*apply(CM1^2, 2, sum))/2+log(initial.y.dist)
      scaled.d <- scale(t(mXtrain) %*% CM1, -intvec, FALSE)
      dimnames(scaled.d) <- list(NULL, names(initial.y.dist))
      Pred.y[i, j, ] <- as.factor(rep(1,nrow(scaled.d)) + scaled.d[,2]>scaled.d[,1])
      
      FinalErr[i, j] <- sum(vytrain != Pred.y[i, j, ])
      scaled.d <- scale(t(mXtest)%*%CM1, -intvec, FALSE)
      dimnames(scaled.d) <- list(NULL, names(initial.y.dist))
      initialXival[i, j, , ] <- scaled.d
      Pred.y.new[i, j, ] <- as.factor(rep(1,nrow(scaled.d)) + scaled.d[,2]>scaled.d[,1])
      
    }
  }
  
  
  
  initialXival <- expit(initialXival[,,,1] - initialXival[,,,2])
  return(initialXival)
  
}

GPDA.sparse.NonStat <- function(mXtrain, vytrain, mXtest, train.cycles = 5, test.cycles = 3, delta){
  
  ntrain <- nrow(mXtrain)
  ntest <- nrow(mXtest)
  n1train <- sum(vytrain)
  n0train <- ntrain - n1train
  p <- ncol(mXtrain)
  
  
  mZ <- matrix(0,nrow=ntrain,ncol=p)
  EDepsInv <- 1/cbind(colVars(mXtrain[vytrain==1,]),colVars(mXtrain[vytrain==0,]),colVars(mXtrain))
  vw <- 1-(0.5+apply(mXtrain,2,function(s){ t.test(s[vytrain==1],s[vytrain==0],var.equal = TRUE)$p.value }))/2
  
  m_tau_nu_inv1 = 1
  m_tau_nu_inv0 = 1
  m_tau_nu_inv = 1
  
  Sigma_Z = array(0,dim=c(p,2,ntrain))
  Sigma_Z[,1,] <- rep(10^(-8),ntrain)
  
  m.tau.inv <- 0.1
  tau2 <- 1
  delta <- 1
  tau_star2 <- 1
  lambda.star <- 6
  mu.lambda <- log(3)
  s2.lambda <- 0.25
  
  lambdaTS <- 1
  alpha <- (log(p)); beta <- 4
  nustar1.est <- c(rep(log(0.5),round(p/2)), rep(log(3),p-round(p/2)-1))
  nustar0.est <- c(rep(log(0.5),round(p/2)), rep(log(3),p-round(p/2)-1))
  nustar.est <- c(rep(log(0.5),round(p/2)), rep(log(3),p-round(p/2)-1))
  SigmaNuStar1.est <- cbind(rep(0.1,p-1),rep(0,p-1))
  SigmaNuStar0.est <- cbind(rep(0.1,p-1),rep(0,p-1))
  SigmaNuStar.est <- cbind(rep(0.1,p-1),rep(0,p-1))
  
  muS = rep(0,p-1)
  zeta = rep(0,ntrain)

  xi <- rep(0.5,ntest)
  zetaPred = rep(0,ntest)
  
  B.nustar <- 4000
  
  for(cyc.count in 1:train.cycles){
    
    if(cyc.count==1){
      
      B.nustar1 = 8000
      
    }
    else{
      
      B.nustar1 = 4000
      
    }
    ###Update mu#######
    MuObj <- FitMu(X=mXtrain,y=vytrain,mZ=mZ,EDepsInv = EDepsInv,vw = vw, NuStarEst1 = nustar1.est, NuStarEst0 = nustar0.est, NuStarEst = nustar.est, VarNuStar1 = SigmaNuStar1.est[,1], VarNuStar0 = SigmaNuStar0.est[,1], VarNuStar = SigmaNuStar.est[,1], m_tau_nu_inv1 = m_tau_nu_inv1, m_tau_nu_inv0 = m_tau_nu_inv0, m_tau_nu_inv = m_tau_nu_inv,delta = delta)
    m_mu1 <- MuObj$m_mu1
    m_mu0 <- MuObj$m_mu0
    m_mu <- MuObj$m_mu
    Sigma_Mu1 <- MuObj$Sigma_Mu1
    Sigma_Mu0 <- MuObj$Sigma_Mu0
    Sigma_Mu <- MuObj$Sigma_Mu

    ###Update tau.star########
    TauStar1Obj <- FitTauStar(m_mu = m_mu1, Sigma_Mu = Sigma_Mu1, nu_est = nustar1.est, VarNuStar = SigmaNuStar1.est[,1], Atau_star = 1,Btau_star = 1,delta = delta)
    r.tau.star1 <- TauStar1Obj[1,1]
    s.tau.star1 <- TauStar1Obj[2,1]
    m_tau_star_inv1 <- r.tau.star1/s.tau.star1
    TauStar0Obj <- FitTauStar(m_mu = m_mu0, Sigma_Mu = Sigma_Mu0, nu_est = nustar0.est, VarNuStar = SigmaNuStar0.est[,1], Atau_star = 1,Btau_star = 1,delta = delta)
    r.tau.star0 <- TauStar0Obj[1,1]
    s.tau.star0 <- TauStar0Obj[2,1]
    m_tau_star_inv0 <- r.tau.star0/s.tau.star0
    TauStarObj <- FitTauStar(m_mu = m_mu, Sigma_Mu = Sigma_Mu, nu_est = nustar.est, VarNuStar = SigmaNuStar.est[,1], Atau_star = 1,Btau_star = 1,delta = delta)
    r.tau.star <- TauStarObj[1,1]
    s.tau.star <- TauStarObj[2,1]
    m_tau_star_inv <- r.tau.star/s.tau.star

    ###Update nu.star#########
    SVBnustar1Obj <- PrecSVBforNuCpp(mZi= m_mu1, SigmaZi = Sigma_Mu1, m_tau_invTS = m_tau_star_inv1, delta=delta, tau2 = tau2, lambda = lambda.star, mStart = rnorm(p-1),B=B.nustar1)
    nustar1.est <- SVBnustar1Obj$mStore
    QNuStar1 <- cbind(c(SVBnustar1Obj$OmegaStore[1,1]^2, SVBnustar1Obj$OmegaStore[1:(p-2),2]^2+ SVBnustar1Obj$OmegaStore[2:(p-1),1]^2) ,c(SVBnustar1Obj$OmegaStore[1:(p-2),1]*SVBnustar1Obj$OmegaStore[1:(p-2),2],0))
    SigmaNuStar1.est <- BandedCholToInvNonStat2(QNuStar1[,1],QNuStar1[,2])
    if(cyc.count==1){
      
      nustar0.est <- nustar1.est
      nustar.est <- nustar1.est
      
    }
    SVBnustar0Obj <- PrecSVBforNuCpp(mZi= m_mu0, SigmaZi = Sigma_Mu0, m_tau_invTS = m_tau_star_inv0, delta=delta, tau2 = tau2, lambda = lambda.star, mStart = nustar0.est,B=B.nustar)
    nustar0.est <- SVBnustar0Obj$mStore
    QNuStar0 <- cbind(c(SVBnustar0Obj$OmegaStore[1,1]^2, SVBnustar0Obj$OmegaStore[1:(p-2),2]^2+ SVBnustar0Obj$OmegaStore[2:(p-1),1]^2) ,c(SVBnustar0Obj$OmegaStore[1:(p-2),1]*SVBnustar0Obj$OmegaStore[1:(p-2),2],0))
    SigmaNuStar0.est <- BandedCholToInvNonStat2(QNuStar0[,1],QNuStar0[,2])
    
    SVBnustarObj <- PrecSVBforNuCpp(mZi= m_mu, SigmaZi = Sigma_Mu, m_tau_invTS = m_tau_star_inv, delta=delta, tau2 = tau2, lambda = lambda.star, mStart = nustar.est,B=B.nustar)
    nustar.est <- SVBnustarObj$mStore
    QNuStar <- cbind(c(SVBnustarObj$OmegaStore[1,1]^2, SVBnustarObj$OmegaStore[1:(p-2),2]^2+ SVBnustarObj$OmegaStore[2:(p-1),1]^2) ,c(SVBnustarObj$OmegaStore[1:(p-2),1]*SVBnustarObj$OmegaStore[1:(p-2),2],0))
    SigmaNuStar.est <- BandedCholToInvNonStat2(QNuStar[,1],QNuStar[,2])

    ###Update lambda.star#########
    NustarsqTotal <- SigmaNuStar1.est[,1] + nustar1.est^2 + SigmaNuStar0.est[,1] + nustar0.est^2 + SigmaNuStar.est[,1] + nustar.est^2
    NuNustarTotal <- SigmaNuStar1.est[1:(p-2),2] + nustar1.est[1:(p-2)]*nustar1.est[2:(p-1)] + SigmaNuStar0.est[1:(p-2),2] + nustar0.est[1:(p-2)]*nustar0.est[2:(p-1)] + SigmaNuStar.est[1:(p-2),2] + nustar.est[1:(p-2)]*nustar.est[2:(p-1)]
    OptLambdaStarObj <- OptimizeLambda(ENusqTotal = NustarsqTotal, ENuNuTotal = NuNustarTotal, m_tau2_inv = (1/tau2), mu.lambda = mu.lambda, s2.lambda = s2.lambda, samp.size = 3, delta = delta)
    lambda.star <- exp(OptLambdaStarObj)

    ###Update epsilon#########
    EpsilonRObj <- FitEpsilonRversion(X=mXtrain,y=vytrain,mZ=mZ, m_mu1 = m_mu1, m_mu0 = m_mu0, m_mu = m_mu, VarMu1 = Sigma_Mu1[,1], VarMu0 = Sigma_Mu0[,1], VarMu = Sigma_Mu[,1], VarZ = t(Sigma_Z[,1,]), vw = vw, Aeps = 0.00001,Beps = 0.000001)
    EDepsInv <- EpsilonRObj$EDepsInv
    ElogSigmaSq <- EpsilonRObj$ELogeps
    #cat("done epsilon","\n")
    if(cyc.count==1){
      
      mS = nustar.est
      Sigma_S = SigmaNuStar.est
      
    }
    
    ###Update Z#########
    ZObj <- FitZ(X = mXtrain, y=vytrain, m_mu1 =m_mu1, m_mu0 = m_mu0, m_mu = m_mu, vw=vw, mEpsilonInvMAT = EDepsInv, zeta = zeta, mS=mS, Sigma_S = Sigma_S, m_tau_inv = m.tau.inv,delta = delta)
    mZ <- ZObj$mZ
    Sigma_Z <- ZObj$SigmaZ
    
    ###Update tauTS#########
    TauTSobj <- FitTauTS2(mZ = mZ, SigmaZ = Sigma_Z, mS = mS, Sigma_S = Sigma_S, zeta=zeta, atauTS = 0.01, btauTS = 0.01,delta=delta)
    m.tau.inv <- TauTSobj$rtauTS/TauTSobj$stauTS
    
    ###Update Gamma#########
    GammaObj <- FitGamma(X=mXtrain, y = vytrain, mZ = mZ, SigmaZ = Sigma_Z, m_mu1 = m_mu1, m_mu0 = m_mu0, m_mu = m_mu, Sigma_Mu1 = Sigma_Mu1, Sigma_Mu0 = Sigma_Mu0, Sigma_Mu = Sigma_Mu , vwOLD = vw, mEpsilonInvMAT = EDepsInv, mLogEpsilonMAT = ElogSigmaSq, delta = 1, alpha = alpha, beta = beta)
    vw <- GammaObj$vw
    
    ###Update alphabeta#########
    sumW <- sum(vw)
    sumAdjW <- sum(vw[1:(p-1)]*vw[2:p])
    
    
    ###Update S###############
    Ezsq <- (mZ^2 + t(Sigma_Z[,1,]))
    Ezz <- (mZ[,1:(p-1)]*mZ[,2:p] + t(Sigma_Z[1:(p-1),2,]))
    ezeta <- exp(zeta)
    repezeta <- rep(ezeta,p)
    
    G1 <- group(c(Ezsq)*repezeta, rep(1:p,rep(ntrain,p)),"sum")
    G2 <- group(c(Ezsq)/repezeta, rep(1:p,rep(ntrain,p)),"sum")
    G3 <- group(c(Ezz)*repezeta[-c(1:ntrain)], rep(1:(p-1),rep(ntrain,p-1)),"sum")
    G4 <- colSums(Ezsq)
    G5 <- colSums(Ezz)
    
    Sobj <- PrecSVBforSCpp(G1=G1, G2=G2, G3=G3, G4=G4, G5=G5, m_tau_inv=m.tau.inv, muS = muS, tau2 = tau2, lambda=lambdaTS, mStart = mS, delta=delta, B= 5100, zeta=zeta)
    mS <- Sobj$mStore
    MainDiagOmegaS <- c(Sobj$OmegaStore[1,1]^2,Sobj$OmegaStore[1:(p-2),2]^2+Sobj$OmegaStore[2:(p-1),1]^2)
    OffDiagOmegaS <- c(Sobj$OmegaStore[1:(p-2),2]*Sobj$OmegaStore[1:(p-2),1],0)
    Sigma_S <- BandedCholToInvNonStat2(MainDiagOmegaS, OffDiagOmegaS)
    
    ###Update lambdaTS###############
    SsqTotal <- Sigma_S[,1] + mS^2
    SSTotal <- Sigma_S[1:(p-2),2] + mS[1:(p-2)]*mS[2:(p-1)]
    OptLambdaTSStarObj <- OptimizeLambda(ENusqTotal = SsqTotal, ENuNuTotal = SSTotal, m_tau2_inv = (1/tau2), mu.lambda = mu.lambda, s2.lambda = s2.lambda, samp.size = 1, delta = 1)
    lambdaTS <- exp(OptLambdaTSStarObj)
    
    ###Update zeta###############
    ZetaObj <- OptimizeZeta3(mZ = mZ, SigmaZ = Sigma_Z, mS = mS, VarS = Sigma_S[,1], sigma2.zeta = 1, m_tau_inv = m.tau.inv, delta = delta)
    zeta = ZetaObj
    
    cat("finished with ", cyc.count, "\n")
  }

  SigmaZPred = array(0,dim=c(p,2,ntest))
  mZPred <- matrix(0,nrow=ntest,ncol=p)
  Sigma_Z1.start <- rowMeans(Sigma_Z[,1,])
  Sigma_Z2.start <- rowMeans(Sigma_Z[,2,])
  for(i in 1:ntest){
    
    SigmaZPred[,1,i] <- Sigma_Z1.start
    SigmaZPred[,2,i] <- Sigma_Z2.start
    
  }
  
  xi <- SetInitialXi(mXtrain = t(mXtrain), vytrain = c(2-vytrain), mXtest = t(mXtest), alpha=0.9, delta=0)
  
  
  for(cyc.count.test in 1:test.cycles){
    
    PredZObj <- FitPredZ(X = mXtest, xi = xi, m_mu1 = m_mu1, m_mu0 = m_mu0, m_mu = m_mu, vw = vw, mEpsilonInvMAT = EDepsInv, zeta = zetaPred, mS=mS, Sigma_S = Sigma_S, m_tau_inv = m.tau.inv, delta = delta)
    mZPred <- PredZObj$mZ
    SigmaZPred <- PredZObj$SigmaZ
    ZetaPredObj <- OptimizeZeta3(mZ = mZPred, SigmaZ = SigmaZPred, mS = mS, VarS = Sigma_S[,1], sigma2.zeta = 1, m_tau_inv = m.tau.inv, delta = delta)
    zetaPred <- ZetaPredObj
    
    xi <- FitPredY(X = mXtest, mZ = mZPred, SigmaZ = SigmaZPred, m_mu1 = m_mu1, m_mu0 = m_mu0, Sigma_Mu1 = Sigma_Mu1, Sigma_Mu0 = Sigma_Mu0, mEpsilonInvMAT = EDepsInv, mLogEpsilonMAT = ElogSigmaSq, vw = vw, delta = delta, logRatio = log(n1train/n0train))
    
  }
  
  return(list(xi=xi,vw=vw,m_mu1 = m_mu1, m_mu0 = m_mu0, m_mu = m_mu, mZ =mZ, Sigma_Z = Sigma_Z, m.tau.inv=m.tau.inv,mS=mS, Sigma_S=Sigma_S, lambdaTS=lambdaTS, zeta=zeta, mZPred=mZPred, SigmaZPred=SigmaZPred, zetaPred = zetaPred))
  
}


Calc.PNmatrix <- function(Predvec, truevec){
  
  if(sum(truevec==2)>0){
    
    truevec <- 2-truevec
    Predvec <- 2-Predvec
    
  }
  TP.id <- which(truevec==1)
  TN.id <- which(truevec==0)
  nFN <- sum(Predvec[TP.id]==0)
  nFP <- sum(Predvec[TN.id]==1)
  nTP <- sum(Predvec==1) - nFP
  nTN <- sum(Predvec==0) - nFN
  
  return(c(nTP,nTN,nFP,nFN))
  
}

Calc.MCC <- function(predvec,truevec){
  
  P.plus <- sum(predvec[which(truevec==1)])
  P.minus <- sum(predvec[which(truevec==0)])
  N.minus <- sum(abs(1-predvec)[which(truevec==0)])
  N.plus <- sum(abs(1-predvec)[which(truevec==1)])
  
  
  MCC.val <- (P.plus*N.minus - P.minus*N.plus)/sqrt((P.plus+P.minus)*(P.plus+N.plus)*(N.minus+P.minus)*(N.minus+N.plus))

  return(MCC.val)
  
}