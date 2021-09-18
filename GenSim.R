
t <- 0:(n.varbs-1)

#Many weak signals
if(sim.set==1){
  
  gamma.true <- rep(0,n.varbs)
  gamma.true[c(51:300,701:950, 1251:1500, 1601:1850, 2001:2250, 3011:3260,3801:4050, 4401:4650)] <- 1
  #ggplot(mapping = aes(x = t, y=gamma.true)) + geom_line()
  sigma.eps.true1 <- gamma.true*exp(sin(0.1*1:n.varbs)) + 0.5*(1-gamma.true)*exp(cos(0.1*1:n.varbs))
  sigma.eps.true0 <- gamma.true*exp(cos(0.1*1:n.varbs)) + 0.5*(1-gamma.true)*exp(cos(0.1*1:n.varbs))
  sigma.eps.true1 <- 0.005*sigma.eps.true1
  sigma.eps.true0 <- 0.005*sigma.eps.true0
  #Please remember to revert
  tau.true <- 4.5; tau2.true <- 2; lambda.true <- 500; tau.star.true <- 2; tau2.star.true <- 2; lambda.star.true <- 500
  #tau.true <- 0.3; tau2.true <- 2; lambda.true <- 500; tau.star.true <- 2; tau2.star.true <- 2; lambda.star.true <- 500
  zeta.train.true <- 0.5*exp((1:nTrain)^0.05) - 1.5; zeta.test.true <- 0.5*exp((1:nTest)^0.05) - 1.5
  LS.star.true <- matrix(0,nrow=(n.varbs-1),ncol=2)
  LS.star.true[,1] <- sqrt(lambda.star.true/(2*delta))
  LS.star.true[1,1] <- 1
  LS.star.true[,2] <- -sqrt(lambda.star.true/(2*delta))*(1-delta/lambda.star.true)
  LS.star.true[n.varbs-1,2] <- 0
  LS.star.true <- 1/sqrt(tau2.star.true)*LS.star.true
  set.seed(10000+sim.set+sim.rep)
  nu.star.true1 <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true1)) + geom_line()
  nu.star.true0 <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true0)) + geom_line()
  nu.star.true <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true)) + geom_line()
  LNS.star.true1 <- matrix(0,nrow=n.varbs,ncol=2)
  LNS.star.true0 <- matrix(0,nrow=n.varbs,ncol=2)
  LNS.star.true <- matrix(0,nrow=n.varbs,ncol=2)
  LNS.star.true1[,1] <- c(1,sqrt(exp(nu.star.true1)/(2*delta)))
  LNS.star.true1[,2] <- c(-sqrt(exp(nu.star.true1)/(2*delta))*(1-delta*exp(-nu.star.true1)), 0)
  LNS.star.true1 <- 1/sqrt(tau.star.true)*LNS.star.true1
  mu.star.true1 <- 0.3*c(ThomasAlgo2(LNS.star.true1[1:(n.varbs-1),2], LNS.star.true1[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
  #ggplot(mapping = aes(x = t, y=mu.star.true1)) + geom_line()
  LNS.star.true0[,1] <- c(1,sqrt(exp(nu.star.true0)/(2*delta)))
  LNS.star.true0[,2] <- c(-sqrt(exp(nu.star.true0)/(2*delta))*(1-delta*exp(-nu.star.true0)), 0)
  LNS.star.true0 <- 1/sqrt(tau.star.true)*LNS.star.true0
  mu.star.true0 <- 0.3*c(ThomasAlgo2(LNS.star.true0[1:(n.varbs-1),2], LNS.star.true0[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
  #ggplot(mapping = aes(x = t, y=mu.star.true0)) + geom_line()
  LNS.star.true[,1] <- c(1,sqrt(exp(nu.star.true)/(2*delta)))
  LNS.star.true[,2] <- c(-sqrt(exp(nu.star.true)/(2*delta))*(1-delta*exp(-nu.star.true)), 0)
  LNS.star.true <- 1/sqrt(tau.star.true)*LNS.star.true
  mu.star.true <- 0.5*c(ThomasAlgo2(LNS.star.true[1:(n.varbs-1),2], LNS.star.true[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
  #ggplot(mapping = aes(x = t, y=mu.star.true)) + geom_line()
  mu.true1 <- (gamma.true*mu.star.true1 + (1-gamma.true)*mu.star.true)
  mu.true0 <- (gamma.true*mu.star.true0 + (1-gamma.true)*mu.star.true)
  #ggplot(mapping = aes(x = t, y=mu.true1)) + geom_line()
  #ggplot(mapping = aes(x = t, y=mu.true0)) + geom_line()
  LS.true <- matrix(0,nrow=(n.varbs-1),ncol=2)
  LS.true[,1] <- sqrt(lambda.true/(2*delta))
  LS.true[1,1] <- 1
  LS.true[,2] <- -sqrt(lambda.true/(2*delta))*(1-delta/lambda.true)
  LS.true[n.varbs-1,2] <- 0
  LS.true <- 1/sqrt(tau2.true)*LS.true
  #S.true <- 7+c(ThomasAlgo2(LS.true[1:(n.varbs-2),2], LS.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  S.true <- 4.5+c(ThomasAlgo2(LS.true[1:(n.varbs-2),2], LS.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  #ggplot(mapping = aes(x = t[2:n.varbs], y=S.true)) + geom_line()
  nu.train.true <- matrix(rep(zeta.train.true,rep(n.varbs-1,nTrain)) + rep(S.true,nTrain),nrow=nTrain,ncol=(n.varbs-1),byrow=T)
  nu.test.true <- matrix(rep(zeta.test.true,rep(n.varbs-1,nTest)) + rep(S.true,nTest),nrow=nTest,ncol=(n.varbs-1),byrow=T)
  mZ.train.true <- matrix(0,nrow=nTrain,ncol=n.varbs)
  #ggplot(mapping = aes(x = t, y=mZ.train.true[1,])) + geom_line()
  mZ.test.true <- matrix(0,nrow=nTest,ncol=n.varbs)
  vy.train <- rbinom(nTrain,1,0.5)
  vy.test <- rbinom(nTest,1,0.5)
  X.train <- matrix(0,nrow=nTrain,ncol=n.varbs)
  X.test <- matrix(0,nrow=nTest,ncol=n.varbs)
  for(i in 1:nTrain){
    
    LNS.true.temp <- matrix(0,nrow=n.varbs,ncol=2)
    LNS.true.temp[,1] <- c(1,sqrt(exp(nu.train.true[i,])/(2*delta)))
    LNS.true.temp[,2] <- c(-sqrt(exp(nu.train.true[i,])/(2*delta))*(1-delta*exp(-nu.train.true[i,])), 0)
    LNS.true.temp <- 1/sqrt(tau.true)*LNS.true.temp
    mZ.train.true[i,] <- c(ThomasAlgo2(LNS.true.temp[1:(n.varbs-1),2], LNS.true.temp[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
    if(vy.train[i]==1){
      
      X.train[i,] <- mZ.train.true[i,] + mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
      
    }
    else{
      
      X.train[i,] <- mZ.train.true[i,] + mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
      
    }
  }
  for(i in 1:nTest){
    
    LNS.true.temp <- matrix(0,nrow=n.varbs,ncol=2)
    LNS.true.temp[,1] <- c(1,sqrt(exp(nu.test.true[i,])/(2*delta)))
    LNS.true.temp[,2] <- c(-sqrt(exp(nu.test.true[i,])/(2*delta))*(1-delta*exp(-nu.test.true[i,])), 0)
    LNS.true.temp <- 1/sqrt(tau.true)*LNS.true.temp
    mZ.test.true[i,] <- c(ThomasAlgo2(LNS.true.temp[1:(n.varbs-1),2], LNS.true.temp[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
    if(vy.test[i]==1){
      
      X.test[i,] <- mZ.test.true[i,] + mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
      
    }
    else{
      
      X.test[i,] <- mZ.test.true[i,] + mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
      
    }
  }
  
}

#Sparse strong signals
if(sim.set==2){
  
  gamma.true <- rep(0,n.varbs)
  gamma.true[c(1601:1725,3201:3375)] <- 1
  #ggplot(mapping = aes(x = t, y=gamma.true)) + geom_line()
  sigma.eps.true1 <- gamma.true*exp(cos(0.1*1:n.varbs)) + 0.5*(1-gamma.true)*exp(cos(0.1*1:n.varbs))
  sigma.eps.true0 <- gamma.true*exp(sin(0.1*1:n.varbs)) + 0.5*(1-gamma.true)*exp(cos(0.1*1:n.varbs))
  sigma.eps.true1 <- 0.5*sigma.eps.true1
  sigma.eps.true0 <- 0.5*sigma.eps.true0
  tau.true <- 1.5; tau2.true <- 2; lambda.true <- 500; tau.star.true <- 2; tau2.star.true <- 2; lambda.star.true <- 500
  zeta.train.true <- 0.5*exp((1:nTrain)^0.05) - 1.5; zeta.test.true <- 0.5*exp((1:nTest)^0.05) - 1.5
  LS.star.true <- matrix(0,nrow=(n.varbs-1),ncol=2)
  LS.star.true[,1] <- sqrt(lambda.star.true/(2*delta))
  LS.star.true[1,1] <- 1
  LS.star.true[,2] <- -sqrt(lambda.star.true/(2*delta))*(1-delta/lambda.star.true)
  LS.star.true[n.varbs-1,2] <- 0
  LS.star.true <- 1/sqrt(tau2.star.true)*LS.star.true
  set.seed(10000+sim.set+sim.rep)
  nu.star.true1 <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true1)) + geom_line()
  nu.star.true0 <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true0)) + geom_line()
  nu.star.true <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true)) + geom_line()
  LNS.star.true1 <- matrix(0,nrow=n.varbs,ncol=2)
  LNS.star.true0 <- matrix(0,nrow=n.varbs,ncol=2)
  LNS.star.true <- matrix(0,nrow=n.varbs,ncol=2)
  LNS.star.true1[,1] <- c(1,sqrt(exp(nu.star.true1)/(2*delta)))
  LNS.star.true1[,2] <- c(-sqrt(exp(nu.star.true1)/(2*delta))*(1-delta*exp(-nu.star.true1)), 0)
  LNS.star.true1 <- 1/sqrt(tau.star.true)*LNS.star.true1
  mu.star.true1 <- 0.7*c(ThomasAlgo2(LNS.star.true1[1:(n.varbs-1),2], LNS.star.true1[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
  #ggplot(mapping = aes(x = t, y=mu.star.true1)) + geom_line()
  LNS.star.true0[,1] <- c(1,sqrt(exp(nu.star.true0)/(2*delta)))
  LNS.star.true0[,2] <- c(-sqrt(exp(nu.star.true0)/(2*delta))*(1-delta*exp(-nu.star.true0)), 0)
  LNS.star.true0 <- 1/sqrt(tau.star.true)*LNS.star.true0
  mu.star.true0 <- 0.7*c(ThomasAlgo2(LNS.star.true0[1:(n.varbs-1),2], LNS.star.true0[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
  #ggplot(mapping = aes(x = t, y=mu.star.true0)) + geom_line()
  LNS.star.true[,1] <- c(1,sqrt(exp(nu.star.true)/(2*delta)))
  LNS.star.true[,2] <- c(-sqrt(exp(nu.star.true)/(2*delta))*(1-delta*exp(-nu.star.true)), 0)
  LNS.star.true <- 1/sqrt(tau.star.true)*LNS.star.true
  mu.star.true <- 0.7*c(ThomasAlgo2(LNS.star.true[1:(n.varbs-1),2], LNS.star.true[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
  #ggplot(mapping = aes(x = t, y=mu.star.true)) + geom_line()
  mu.true1 <- (gamma.true*mu.star.true1 + (1-gamma.true)*mu.star.true)
  mu.true0 <- (gamma.true*mu.star.true0 + (1-gamma.true)*mu.star.true)
  #ggplot(mapping = aes(x = t, y=mu.true1)) + geom_line()
  #ggplot(mapping = aes(x = t, y=mu.true0)) + geom_line()
  LS.true <- matrix(0,nrow=(n.varbs-1),ncol=2)
  LS.true[,1] <- sqrt(lambda.true/(2*delta))
  LS.true[1,1] <- 1
  LS.true[,2] <- -sqrt(lambda.true/(2*delta))*(1-delta/lambda.true)
  LS.true[n.varbs-1,2] <- 0
  LS.true <- 1/sqrt(tau2.true)*LS.true
  #S.true <- 7+c(ThomasAlgo2(LS.true[1:(n.varbs-2),2], LS.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  S.true <- 6+c(ThomasAlgo2(LS.true[1:(n.varbs-2),2], LS.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
  #ggplot(mapping = aes(x = t[2:n.varbs], y=S.true)) + geom_line()
  nu.train.true <- matrix(rep(zeta.train.true,rep(n.varbs-1,nTrain)) + rep(S.true,nTrain),nrow=nTrain,ncol=(n.varbs-1),byrow=T)
  nu.test.true <- matrix(rep(zeta.test.true,rep(n.varbs-1,nTest)) + rep(S.true,nTest),nrow=nTest,ncol=(n.varbs-1),byrow=T)
  mZ.train.true <- matrix(0,nrow=nTrain,ncol=n.varbs)
  #ggplot(mapping = aes(x = t, y=mZ.train.true[1,])) + geom_line()
  mZ.test.true <- matrix(0,nrow=nTest,ncol=n.varbs)
  vy.train <- rbinom(nTrain,1,0.5)
  vy.test <- rbinom(nTest,1,0.5)
  X.train <- matrix(0,nrow=nTrain,ncol=n.varbs)
  X.test <- matrix(0,nrow=nTest,ncol=n.varbs)
  for(i in 1:nTrain){
    
    LNS.true.temp <- matrix(0,nrow=n.varbs,ncol=2)
    LNS.true.temp[,1] <- c(1,sqrt(exp(nu.train.true[i,])/(2*delta)))
    LNS.true.temp[,2] <- c(-sqrt(exp(nu.train.true[i,])/(2*delta))*(1-delta*exp(-nu.train.true[i,])), 0)
    LNS.true.temp <- 1/sqrt(tau.true)*LNS.true.temp
    mZ.train.true[i,] <- c(ThomasAlgo2(LNS.true.temp[1:(n.varbs-1),2], LNS.true.temp[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
    if(vy.train[i]==1){
      
      X.train[i,] <- mZ.train.true[i,] + mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
      
    }
    else{
      
      X.train[i,] <- mZ.train.true[i,] + mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
      
    }
  }
  for(i in 1:nTest){
    
    LNS.true.temp <- matrix(0,nrow=n.varbs,ncol=2)
    LNS.true.temp[,1] <- c(1,sqrt(exp(nu.test.true[i,])/(2*delta)))
    LNS.true.temp[,2] <- c(-sqrt(exp(nu.test.true[i,])/(2*delta))*(1-delta*exp(-nu.test.true[i,])), 0)
    LNS.true.temp <- 1/sqrt(tau.true)*LNS.true.temp
    mZ.test.true[i,] <- c(ThomasAlgo2(LNS.true.temp[1:(n.varbs-1),2], LNS.true.temp[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
    if(vy.test[i]==1){
      
      X.test[i,] <- mZ.test.true[i,] + mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
      
    }
    else{
      
      X.test[i,] <- mZ.test.true[i,] + mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
      
    }
  }
  
}


#Independence
if(sim.set==3){
  
  gamma.true <- rep(0,n.varbs)
  gamma.true[rep(seq(201,450,1),2) + rep(c(1000,2000),c(250,250))] <- 1
  set.seed(987654321)
  mu.true.star1 <- 0.25*rt(n.varbs,df=5)
  mu.true.star0 <- 0.25*rt(n.varbs,df=5)
  mu.true.star <- rt(n.varbs,df=5)
 
  
  mu.true1 <- gamma.true*mu.true.star1 + (1-gamma.true)*mu.true.star
  mu.true0 <- gamma.true*mu.true.star0 + (1-gamma.true)*mu.true.star
  
  sigma.eps.true1 <- 2 + cos(0.01*1:n.varbs)
  sigma.eps.true0 <- sigma.eps.true1
  set.seed(sim.rep + sim.set)
  vy.train <- rbinom(nTrain,1,0.5)
  vy.test <- rbinom(nTest,1,0.5)
  X.train <- matrix(0,nrow=nTrain,ncol=n.varbs)
  X.test <- matrix(0,nrow=nTest,ncol=n.varbs)
  for(i in 1:nTrain){
    
    if(vy.train[i]==1){
      
      X.train[i,] <- mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
      
    }
    else{
      
      X.train[i,] <- mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
      
    }
  }
  
  for(i in 1:nTest){
    
    if(vy.test[i]==1){
      
      X.test[i,] <- mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
      
    }
    else{
      
      X.test[i,] <- mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
      
    }
  }
  
}

#Uniform covariance structure
if(sim.set==4){
  
  gamma.true <- rep(0,n.varbs)
  gamma.true[rep(seq(401,700,1),2) + rep(c(1000,3000),c(300,300))] <- 1
  mu.true.star1 <- 0.085*(1:n.varbs)^0.1
  mu.true.star0 <- -0.085*(1:n.varbs)^0.2
  mu.true.star <- sin(0.005*1:n.varbs)
  
  mu.true1 <- gamma.true*mu.true.star1 + (1-gamma.true)*mu.true.star
  mu.true0 <- gamma.true*mu.true.star0 + (1-gamma.true)*mu.true.star
  
  set.seed(10000+sim.rep + sim.set)
  vy.train <- rbinom(nTrain,1,0.5)
  vy.test <- sort(rbinom(nTest,1,0.5))
  X.train <- matrix(0,nrow=nTrain,ncol=n.varbs)
  #X.test <- matrix(0,nrow=nTest,ncol=n.varbs)
  
  for(i in 1:nTrain){
    
    if(vy.train[i]==1){
      
      X.train[i,] <- mu.true1 + cholM %*% rnorm(n.varbs)
      
    }
    else{
      
      X.train[i,] <- mu.true0 + cholM %*% rnorm(n.varbs)
      
    }
  }
  
  num1 <- sum(vy.test)
  num0 <- sum(vy.test==0)
  X.test0 <- matrix(rep(mu.true0,num0) +  c(cholM %*% matrix(rnorm(n.varbs*num0),nrow=n.varbs,ncol=num0)), nrow=n.varbs,ncol=num0)
  X.test1 <- matrix(rep(mu.true1,num1) +  c(cholM %*% matrix(rnorm(n.varbs*num1),nrow=n.varbs,ncol=num1)), nrow=n.varbs,ncol=num1)
  X.test <- rbind(t(X.test0), t(X.test1))
  
    
}





# #Many weak signals
# if(sim.set==1){
#   
#   gamma.true <- rep(0,n.varbs)
#   gamma.true[c(51:300,701:950, 1251:1500, 1601:1850, 2001:2250, 3011:3260,3801:4050, 4401:4650)] <- 1
#   #ggplot(mapping = aes(x = t, y=gamma.true)) + geom_line()
#   sigma.eps.true1 <- gamma.true*exp(sin(0.1*1:n.varbs)) + 0.5*(1-gamma.true)*exp(cos(0.1*1:n.varbs))
#   sigma.eps.true0 <- gamma.true*exp(cos(0.1*1:n.varbs)) + 0.5*(1-gamma.true)*exp(cos(0.1*1:n.varbs))
#   sigma.eps.true1 <- 0.005*sigma.eps.true1
#   sigma.eps.true0 <- 0.005*sigma.eps.true0
#   #Please remember to revert
#   tau.true <- 4.5; tau2.true <- 2; lambda.true <- 500; tau.star.true <- 2; tau2.star.true <- 2; lambda.star.true <- 500
#   #tau.true <- 0.3; tau2.true <- 2; lambda.true <- 500; tau.star.true <- 2; tau2.star.true <- 2; lambda.star.true <- 500
#   zeta.train.true <- 0.5*exp((1:nTrain)^0.05) - 1.5; zeta.test.true <- 0.5*exp((1:nTest)^0.05) - 1.5
#   LS.star.true <- matrix(0,nrow=(n.varbs-1),ncol=2)
#   LS.star.true[,1] <- sqrt(lambda.star.true/(2*delta))
#   LS.star.true[1,1] <- 1
#   LS.star.true[,2] <- -sqrt(lambda.star.true/(2*delta))*(1-delta/lambda.star.true)
#   LS.star.true[n.varbs-1,2] <- 0
#   LS.star.true <- 1/sqrt(tau2.star.true)*LS.star.true
#   set.seed(10000+sim.set+sim.rep)
#   nu.star.true1 <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true1)) + geom_line()
#   nu.star.true0 <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true0)) + geom_line()
#   nu.star.true <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true)) + geom_line()
#   LNS.star.true1 <- matrix(0,nrow=n.varbs,ncol=2)
#   LNS.star.true0 <- matrix(0,nrow=n.varbs,ncol=2)
#   LNS.star.true <- matrix(0,nrow=n.varbs,ncol=2)
#   LNS.star.true1[,1] <- c(1,sqrt(exp(nu.star.true1)/(2*delta)))
#   LNS.star.true1[,2] <- c(-sqrt(exp(nu.star.true1)/(2*delta))*(1-delta*exp(-nu.star.true1)), 0)
#   LNS.star.true1 <- 1/sqrt(tau.star.true)*LNS.star.true1
#   mu.star.true1 <- 0.3*c(ThomasAlgo2(LNS.star.true1[1:(n.varbs-1),2], LNS.star.true1[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#   #ggplot(mapping = aes(x = t, y=mu.star.true1)) + geom_line()
#   LNS.star.true0[,1] <- c(1,sqrt(exp(nu.star.true0)/(2*delta)))
#   LNS.star.true0[,2] <- c(-sqrt(exp(nu.star.true0)/(2*delta))*(1-delta*exp(-nu.star.true0)), 0)
#   LNS.star.true0 <- 1/sqrt(tau.star.true)*LNS.star.true0
#   mu.star.true0 <- 0.3*c(ThomasAlgo2(LNS.star.true0[1:(n.varbs-1),2], LNS.star.true0[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#   #ggplot(mapping = aes(x = t, y=mu.star.true0)) + geom_line()
#   LNS.star.true[,1] <- c(1,sqrt(exp(nu.star.true)/(2*delta)))
#   LNS.star.true[,2] <- c(-sqrt(exp(nu.star.true)/(2*delta))*(1-delta*exp(-nu.star.true)), 0)
#   LNS.star.true <- 1/sqrt(tau.star.true)*LNS.star.true
#   mu.star.true <- 0.5*c(ThomasAlgo2(LNS.star.true[1:(n.varbs-1),2], LNS.star.true[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#   #ggplot(mapping = aes(x = t, y=mu.star.true)) + geom_line()
#   mu.true1 <- (gamma.true*mu.star.true1 + (1-gamma.true)*mu.star.true)
#   mu.true0 <- (gamma.true*mu.star.true0 + (1-gamma.true)*mu.star.true)
#   #ggplot(mapping = aes(x = t, y=mu.true1)) + geom_line()
#   #ggplot(mapping = aes(x = t, y=mu.true0)) + geom_line()
#   LS.true <- matrix(0,nrow=(n.varbs-1),ncol=2)
#   LS.true[,1] <- sqrt(lambda.true/(2*delta))
#   LS.true[1,1] <- 1
#   LS.true[,2] <- -sqrt(lambda.true/(2*delta))*(1-delta/lambda.true)
#   LS.true[n.varbs-1,2] <- 0
#   LS.true <- 1/sqrt(tau2.true)*LS.true
#   #S.true <- 7+c(ThomasAlgo2(LS.true[1:(n.varbs-2),2], LS.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   S.true <- 4.5+c(ThomasAlgo2(LS.true[1:(n.varbs-2),2], LS.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   #ggplot(mapping = aes(x = t[2:n.varbs], y=S.true)) + geom_line()
#   nu.train.true <- matrix(rep(zeta.train.true,rep(n.varbs-1,nTrain)) + rep(S.true,nTrain),nrow=nTrain,ncol=(n.varbs-1),byrow=T)
#   nu.test.true <- matrix(rep(zeta.test.true,rep(n.varbs-1,nTest)) + rep(S.true,nTest),nrow=nTest,ncol=(n.varbs-1),byrow=T)
#   mZ.train.true <- matrix(0,nrow=nTrain,ncol=n.varbs)
#   #ggplot(mapping = aes(x = t, y=mZ.train.true[1,])) + geom_line()
#   mZ.test.true <- matrix(0,nrow=nTest,ncol=n.varbs)
#   vy.train <- rbinom(nTrain,1,0.5)
#   vy.test <- rbinom(nTest,1,0.5)
#   X.train <- matrix(0,nrow=nTrain,ncol=n.varbs)
#   X.test <- matrix(0,nrow=nTest,ncol=n.varbs)
#   for(i in 1:nTrain){
#     
#     LNS.true.temp <- matrix(0,nrow=n.varbs,ncol=2)
#     LNS.true.temp[,1] <- c(1,sqrt(exp(nu.train.true[i,])/(2*delta)))
#     LNS.true.temp[,2] <- c(-sqrt(exp(nu.train.true[i,])/(2*delta))*(1-delta*exp(-nu.train.true[i,])), 0)
#     LNS.true.temp <- 1/sqrt(tau.true)*LNS.true.temp
#     mZ.train.true[i,] <- c(ThomasAlgo2(LNS.true.temp[1:(n.varbs-1),2], LNS.true.temp[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#     if(vy.train[i]==1){
#       
#       X.train[i,] <- mZ.train.true[i,] + mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
#       
#     }
#     else{
#       
#       X.train[i,] <- mZ.train.true[i,] + mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
#       
#     }
#   }
#   for(i in 1:nTest){
#     
#     LNS.true.temp <- matrix(0,nrow=n.varbs,ncol=2)
#     LNS.true.temp[,1] <- c(1,sqrt(exp(nu.test.true[i,])/(2*delta)))
#     LNS.true.temp[,2] <- c(-sqrt(exp(nu.test.true[i,])/(2*delta))*(1-delta*exp(-nu.test.true[i,])), 0)
#     LNS.true.temp <- 1/sqrt(tau.true)*LNS.true.temp
#     mZ.test.true[i,] <- c(ThomasAlgo2(LNS.true.temp[1:(n.varbs-1),2], LNS.true.temp[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#     if(vy.test[i]==1){
#       
#       X.test[i,] <- mZ.test.true[i,] + mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
#       
#     }
#     else{
#       
#       X.test[i,] <- mZ.test.true[i,] + mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
#       
#     }
#   }
#   
# }
# 
# 
# 
# #Many weak signals
# if(sim.set==7){
#   
#   gamma.true <- rep(0,n.varbs)
#   gamma.true[c(51:300,701:950, 1251:1500, 1601:1850, 2001:2250, 3011:3260,3801:4050, 4401:4650)] <- 1
#   #ggplot(mapping = aes(x = t, y=gamma.true)) + geom_line()
#   sigma.eps.true1 <- gamma.true*exp(sin(0.1*1:n.varbs)) + 0.5*(1-gamma.true)*exp(cos(0.1*1:n.varbs))
#   sigma.eps.true0 <- gamma.true*exp(cos(0.1*1:n.varbs)) + 0.5*(1-gamma.true)*exp(cos(0.1*1:n.varbs))
#   sigma.eps.true1 <- 0.005*sigma.eps.true1
#   sigma.eps.true0 <- 0.005*sigma.eps.true0
#   #Please remember to revert
#   tau.true <- 1; tau2.true <- 2; lambda.true <- 200; tau.star.true <- 2; tau2.star.true <- 2; lambda.star.true <- 500
#   #tau.true <- 0.3; tau2.true <- 2; lambda.true <- 500; tau.star.true <- 2; tau2.star.true <- 2; lambda.star.true <- 500
#   zeta.train.true <- 0.5*exp((1:nTrain)^0.05) - 1.5; zeta.test.true <- 0.5*exp((1:nTest)^0.05) - 1.5
#   LS.star.true <- matrix(0,nrow=(n.varbs-1),ncol=2)
#   LS.star.true[,1] <- sqrt(lambda.star.true/(2*delta))
#   LS.star.true[1,1] <- 1
#   LS.star.true[,2] <- -sqrt(lambda.star.true/(2*delta))*(1-delta/lambda.star.true)
#   LS.star.true[n.varbs-1,2] <- 0
#   LS.star.true <- 1/sqrt(tau2.star.true)*LS.star.true
#   set.seed(10000+sim.set+sim.rep)
#   nu.star.true1 <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true1)) + geom_line()
#   nu.star.true0 <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true0)) + geom_line()
#   nu.star.true <- 4+c(ThomasAlgo2(LS.star.true[1:(n.varbs-2),2], LS.star.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   #ggplot(mapping = aes(x = t[2:n.varbs], y=nu.star.true)) + geom_line()
#   LNS.star.true1 <- matrix(0,nrow=n.varbs,ncol=2)
#   LNS.star.true0 <- matrix(0,nrow=n.varbs,ncol=2)
#   LNS.star.true <- matrix(0,nrow=n.varbs,ncol=2)
#   LNS.star.true1[,1] <- c(1,sqrt(exp(nu.star.true1)/(2*delta)))
#   LNS.star.true1[,2] <- c(-sqrt(exp(nu.star.true1)/(2*delta))*(1-delta*exp(-nu.star.true1)), 0)
#   LNS.star.true1 <- 1/sqrt(tau.star.true)*LNS.star.true1
#   mu.star.true1 <- 0.3*c(ThomasAlgo2(LNS.star.true1[1:(n.varbs-1),2], LNS.star.true1[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#   #ggplot(mapping = aes(x = t, y=mu.star.true1)) + geom_line()
#   LNS.star.true0[,1] <- c(1,sqrt(exp(nu.star.true0)/(2*delta)))
#   LNS.star.true0[,2] <- c(-sqrt(exp(nu.star.true0)/(2*delta))*(1-delta*exp(-nu.star.true0)), 0)
#   LNS.star.true0 <- 1/sqrt(tau.star.true)*LNS.star.true0
#   mu.star.true0 <- 0.3*c(ThomasAlgo2(LNS.star.true0[1:(n.varbs-1),2], LNS.star.true0[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#   #ggplot(mapping = aes(x = t, y=mu.star.true0)) + geom_line()
#   LNS.star.true[,1] <- c(1,sqrt(exp(nu.star.true)/(2*delta)))
#   LNS.star.true[,2] <- c(-sqrt(exp(nu.star.true)/(2*delta))*(1-delta*exp(-nu.star.true)), 0)
#   LNS.star.true <- 1/sqrt(tau.star.true)*LNS.star.true
#   mu.star.true <- 0.5*c(ThomasAlgo2(LNS.star.true[1:(n.varbs-1),2], LNS.star.true[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#   #ggplot(mapping = aes(x = t, y=mu.star.true)) + geom_line()
#   mu.true1 <- (gamma.true*mu.star.true1 + (1-gamma.true)*mu.star.true)
#   mu.true0 <- (gamma.true*mu.star.true0 + (1-gamma.true)*mu.star.true)
#   #ggplot(mapping = aes(x = t, y=mu.true1)) + geom_line()
#   #ggplot(mapping = aes(x = t, y=mu.true0)) + geom_line()
#   LS.true <- matrix(0,nrow=(n.varbs-1),ncol=2)
#   LS.true[,1] <- sqrt(lambda.true/(2*delta))
#   LS.true[1,1] <- 1
#   LS.true[,2] <- -sqrt(lambda.true/(2*delta))*(1-delta/lambda.true)
#   LS.true[n.varbs-1,2] <- 0
#   LS.true <- 1/sqrt(tau2.true)*LS.true
#   #S.true <- 7+c(ThomasAlgo2(LS.true[1:(n.varbs-2),2], LS.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   S.true <- 4.5+c(ThomasAlgo2(LS.true[1:(n.varbs-2),2], LS.true[,1], rep(0,n.varbs-2), rnorm(n.varbs-1)))
#   #ggplot(mapping = aes(x = t[2:n.varbs], y=S.true)) + geom_line()
#   nu.train.true <- matrix(rep(zeta.train.true,rep(n.varbs-1,nTrain)) + rep(S.true,nTrain),nrow=nTrain,ncol=(n.varbs-1),byrow=T)
#   nu.test.true <- matrix(rep(zeta.test.true,rep(n.varbs-1,nTest)) + rep(S.true,nTest),nrow=nTest,ncol=(n.varbs-1),byrow=T)
#   mZ.train.true <- matrix(0,nrow=nTrain,ncol=n.varbs)
#   #ggplot(mapping = aes(x = t, y=mZ.train.true[1,])) + geom_line()
#   mZ.test.true <- matrix(0,nrow=nTest,ncol=n.varbs)
#   vy.train <- rbinom(nTrain,1,0.5)
#   vy.test <- rbinom(nTest,1,0.5)
#   X.train <- matrix(0,nrow=nTrain,ncol=n.varbs)
#   X.test <- matrix(0,nrow=nTest,ncol=n.varbs)
#   for(i in 1:nTrain){
#     
#     LNS.true.temp <- matrix(0,nrow=n.varbs,ncol=2)
#     LNS.true.temp[,1] <- c(1,sqrt(exp(nu.train.true[i,])/(2*delta)))
#     LNS.true.temp[,2] <- c(-sqrt(exp(nu.train.true[i,])/(2*delta))*(1-delta*exp(-nu.train.true[i,])), 0)
#     LNS.true.temp <- 1/sqrt(tau.true)*LNS.true.temp
#     mZ.train.true[i,] <- c(ThomasAlgo2(LNS.true.temp[1:(n.varbs-1),2], LNS.true.temp[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#     if(vy.train[i]==1){
#       
#       X.train[i,] <- mZ.train.true[i,] + mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
#       
#     }
#     else{
#       
#       X.train[i,] <- mZ.train.true[i,] + mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
#       
#     }
#   }
#   for(i in 1:nTest){
#     
#     LNS.true.temp <- matrix(0,nrow=n.varbs,ncol=2)
#     LNS.true.temp[,1] <- c(1,sqrt(exp(nu.test.true[i,])/(2*delta)))
#     LNS.true.temp[,2] <- c(-sqrt(exp(nu.test.true[i,])/(2*delta))*(1-delta*exp(-nu.test.true[i,])), 0)
#     LNS.true.temp <- 1/sqrt(tau.true)*LNS.true.temp
#     mZ.test.true[i,] <- c(ThomasAlgo2(LNS.true.temp[1:(n.varbs-1),2], LNS.true.temp[,1], rep(0,n.varbs-1), rnorm(n.varbs)))
#     if(vy.test[i]==1){
#       
#       X.test[i,] <- mZ.test.true[i,] + mu.true1 + rnorm(n.varbs)*sqrt(sigma.eps.true1)
#       
#     }
#     else{
#       
#       X.test[i,] <- mZ.test.true[i,] + mu.true0 + rnorm(n.varbs)*sqrt(sigma.eps.true0)
#       
#     }
#   }
#   
# }