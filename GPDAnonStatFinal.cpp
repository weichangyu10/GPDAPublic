#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::vec expit(const arma::vec& v){
  
  arma::vec ans(v.n_elem);
  for(int i = 0; i < v.n_elem; i++){
    
    ans(i) = 1/(1+exp(-v(i)));
    
  }
  return ans;
  
}

// [[Rcpp::export]]
arma::mat EfficientCholNonStat2(const arma::vec& maindiag, const arma::vec& offdiag){
  
  const int T = maindiag.n_elem;
  arma::mat L(T,2);
  L.zeros();
  double s;
  int j;
  
  s = pow(L(0,0),2.0);
  L(0,0) = sqrt(maindiag(0) - s); 
  
  for(int i = 1; i < T; i++){
    
    j = i-1;
    s = L(i,1)*L(i-1,0);
    L(i,1) = (offdiag(j) - s)/L(j,0);
    
    j = i;
    s = L(i,1)*L(i,1);
    L(i,0) = sqrt(maindiag(i) - s);
    
    L(i-1,1) = L(i,1);
    
  }
  
  L(T-1,1) = 0.0;
  return L;
  
  
}

// [[Rcpp::export]]
arma::mat BandedCholToInvNonStat2(const arma::vec& maindiag, const arma::vec& offdiag){
  
  const int T = maindiag.n_elem;
  arma::mat Sym(T,3);
  Sym.zeros();
  int i;
  
  arma::mat L = EfficientCholNonStat2(maindiag, offdiag);
  Sym(T-1,1) += 1/pow(L(T-1,0),2.0);
  Sym(T-2,2) = -L(T-2,1)*Sym(T-1,1)/L(T-2,0);
  Sym(T-1,0) = Sym(T-2,2);
  for(int j = (T-2); j >= 1; j--){
    
    i = j;
    Sym(i,1) = - L(i,1)*Sym(j+1,0)/L(i,0);
    Sym(i,1) += 1/pow(L(i,0),2);
    
    i = j - 1;
    Sym(i,2) = -L(i,1)*Sym(j,1)/L(i,0);
    Sym(j,0) = Sym(i,2);


    
  }

  Sym(0,1) = - L(0,1)*Sym(1,0)/L(0,0);
  Sym(0,1) += 1/pow(L(0,0),2);
  
  arma:: mat ansMAT(T,2);
  ansMAT = Sym.cols(1,2);
  //col 1 (length T) of ansMAT is main diagonal of the solution. Col 2 (length T-1) is first off-diagonal.
  return ansMAT;
  
}

// [[Rcpp::export]]
arma::vec ThomasAlgo2(const arma::vec& va, const arma::vec& vb, const arma::vec& vc, const arma::vec& vd){
  
  const int n = vd.n_elem;
  arma::vec cprime(n-1);
  arma::vec dprime(n);
  arma::vec vx(n);
  //va is lower off-diag
  //vb is main diag
  //vc is upper off-diag
  //vd is RHS
  cprime(0) = vc(0)/vb(0);
  dprime(0) = vd(0)/vb(0);
  for(int i = 1; i < (n-1) ; i++){
    
    cprime(i) = vc(i)/(vb(i) - va(i-1)*cprime(i-1));
    dprime(i) = (vd(i) - va(i-1)*dprime(i-1))/(vb(i) - va(i-1)*cprime(i-1));
    
  }
  dprime(n-1) = (vd(n-1) - va(n-2)*dprime(n-2))/(vb(n-1) - va(n-2)*cprime(n-2));
  vx(n-1) = dprime(n-1);
  for(int j = (n-2); j >= 0; j--){
    
    vx(j) = dprime(j) - cprime(j)*vx(j+1);
    
  }
  
  return vx;
}

//[[Rcpp::export]]
double CalculateK(const arma::vec& m_mu1, const arma::vec& m_mu0, const arma::vec& m_mu, const arma::mat& Sigma_Mu1, const arma::mat& Sigma_Mu0, const arma::mat& Sigma_Mu, const double lambda_star, const double& delta){
  
  const int p = m_mu1.n_elem;
  const double a0 = 1/sqrt(8)*( sqrt(delta/lambda_star) + sqrt(delta/lambda_star + 4*lambda_star/delta));
  const double a1 = 1/sqrt(8)*( sqrt(delta/lambda_star) - sqrt(delta/lambda_star + 4*lambda_star/delta));
  double S1 = 0;
  double S2 = 0;
  double S3 = 0;
  for(int t = 0; t < p; t++){
    
    if(t==0){
      
      S1 += pow(m_mu1(t),2) + pow(m_mu0(t),2) + pow(m_mu(t),2) + Sigma_Mu1(t,0) + Sigma_Mu0(t,0) + Sigma_Mu(t,0);
      S3 += m_mu1(t)*m_mu1(t+1) + m_mu0(t)*m_mu0(t+1) + m_mu(t)*m_mu(t+1) + Sigma_Mu1(t,1) + Sigma_Mu0(t,1) + Sigma_Mu(t,1);
      
    }
    else{
      
      if(t==(p-1)){
        
        S2 += pow(m_mu1(t),2) + pow(m_mu0(t),2) + pow(m_mu(t),2) + Sigma_Mu1(t,0) + Sigma_Mu0(t,0) + Sigma_Mu(t,0);
        
      }
      else{
        
        S1 += pow(m_mu1(t),2) + pow(m_mu0(t),2) + pow(m_mu(t),2) + Sigma_Mu1(t,0) + Sigma_Mu0(t,0) + Sigma_Mu(t,0);
        S2 += pow(m_mu1(t),2) + pow(m_mu0(t),2) + pow(m_mu(t),2) + Sigma_Mu1(t,0) + Sigma_Mu0(t,0) + Sigma_Mu(t,0);
        S3 += m_mu1(t)*m_mu1(t+1) + m_mu0(t)*m_mu0(t+1) + m_mu(t)*m_mu(t+1) + Sigma_Mu1(t,1) + Sigma_Mu0(t,1) + Sigma_Mu(t,1);
        
      }
      
      
      
    }
   
    
    
  }
  double ans = pow(a0, 2)*S1 + pow(a1, 2)*S2 + 2*a0*a1*S3;
  return ans;
  
}


// [[Rcpp::export]]
List FitMu(const arma::mat& X, const arma::vec& y, const arma::mat& mZ, const arma::mat& EDepsInv, const arma::vec& vw, const arma::vec& NuStarEst1, const arma::vec& NuStarEst0, const arma::vec& NuStarEst, const arma::vec& VarNuStar1, const arma::vec& VarNuStar0, const arma::vec& VarNuStar, const double m_tau_nu_inv1, const double m_tau_nu_inv0, const double m_tau_nu_inv, const double& delta){

  const int p = X.n_cols;
  const int n = y.n_elem;
  const int n1 = accu(y);
  const int n0 = n - n1;
  const arma::vec onesVec = ones(p);
  const arma::vec halfvec = 0.5*onesVec.rows(0,p-2);
  
  const arma::vec Eb0_sq1 = 0.5*exp(NuStarEst1 + 0.5*VarNuStar1)/delta;
  const arma::vec Eb1_sq1 = Eb0_sq1 + 0.5*delta*exp(-NuStarEst1 + 0.5*VarNuStar1) - onesVec.rows(0,p-2);
  const arma::vec Eb0b1_1 = halfvec - Eb0_sq1;
  
  const arma::vec Eb0_sq0 = 0.5*exp(NuStarEst0 + 0.5*VarNuStar0)/delta;
  const arma::vec Eb1_sq0 = Eb0_sq0 + 0.5*delta*exp(-NuStarEst0 + 0.5*VarNuStar0) - onesVec.rows(0,p-2);
  const arma::vec Eb0b1_0 = halfvec - Eb0_sq0;
  
  const arma::vec Eb0_sq = 0.5*exp(NuStarEst + 0.5*VarNuStar)/delta;
  const arma::vec Eb1_sq = Eb0_sq + 0.5*delta*exp(-NuStarEst + 0.5*VarNuStar) - onesVec.rows(0,p-2);
  const arma::vec Eb0b1 = halfvec - Eb0_sq;
  
  //col 0 = lower off-diag, col 1 = main diag
  arma::mat Qmu1(p,2);
  Qmu1.zeros();
  arma::mat Qmu0(p,2);
  Qmu0.zeros();
  arma::mat Qmu(p,2);
  Qmu.zeros();
  Qmu1.submat(0,0,p-2,0) =  Eb0b1_1;
  Qmu1.submat(0,1,p-2,1) = Eb1_sq1;
  Qmu1.submat(1,1,p-1,1) = Qmu1.submat(1,1,p-1,1) + Eb0_sq1;
  Qmu1(0,1) = Qmu1(0,1) + 1.0;
  Qmu1 = m_tau_nu_inv1*Qmu1;
  //EDeps col 0: k=1
  Qmu1.col(1) = n1*vw%EDepsInv.col(0) + Qmu1.col(1);
  Qmu0.submat(0,0,p-2,0) =  Eb0b1_0;
  Qmu0.submat(0,1,p-2,1) = Eb1_sq0;
  Qmu0.submat(1,1,p-1,1) = Qmu0.submat(1,1,p-1,1) + Eb0_sq0;
  Qmu0(0,1) = Qmu0(0,1) + 1.0;
  Qmu0 = m_tau_nu_inv0*Qmu0;
  //EDeps col 1: k=0
  Qmu0.col(1) = n0*vw%EDepsInv.col(1) + Qmu0.col(1);
  Qmu.submat(0,0,p-2,0) =  Eb0b1;
  Qmu.submat(0,1,p-2,1) = Eb1_sq;
  Qmu.submat(1,1,p-1,1) = Qmu.submat(1,1,p-1,1) + Eb0_sq;
  Qmu(0,1) = Qmu(0,1) + 1.0;
  Qmu = m_tau_nu_inv*Qmu;
  Qmu.col(1) = n*(ones(p)-vw)%EDepsInv.col(2) + Qmu.col(1);
  const uvec inds1 = find(y == 1);
  const uvec inds0 = find(y == 0);
  const arma::mat Xtrain1 =  X.rows(inds1);
  const arma::mat Xtrain0 =  X.rows(inds0);
  const arma::mat mZ1 =  mZ.rows(inds1);
  const arma::mat mZ0 =  mZ.rows(inds0);
  
  arma::vec wtempD1 = EDepsInv.col(0)%(vw%(sum( (Xtrain1 - mZ1), 0)).t());
  arma::vec m_mu1 = ThomasAlgo2(Qmu1.submat(0,0,p-2,0), Qmu1.col(1), Qmu1.submat(0,0,p-2,0), wtempD1);

  arma::vec wtempD0 = EDepsInv.col(1)%(vw%(sum( (Xtrain0 - mZ0), 0)).t());
  arma::vec m_mu0 = ThomasAlgo2(Qmu0.submat(0,0,p-2,0), Qmu0.col(1), Qmu0.submat(0,0,p-2,0), wtempD0);
  
  arma::vec wtempD = EDepsInv.col(2)%( (ones(p)-vw)%(sum( (X - mZ), 0)).t());
  arma::vec m_mu = ThomasAlgo2(Qmu.submat(0,0,p-2,0), Qmu.col(1), Qmu.submat(0,0,p-2,0), wtempD);
  
  arma::mat Sigma_Mu1 = BandedCholToInvNonStat2(Qmu1.col(1), Qmu1.col(0));
  arma::mat Sigma_Mu0 = BandedCholToInvNonStat2(Qmu0.col(1), Qmu0.col(0));
  arma::mat Sigma_Mu = BandedCholToInvNonStat2(Qmu.col(1), Qmu.col(0));
  
  return Rcpp::List::create(Rcpp::Named("m_mu1") = m_mu1, Rcpp::Named("m_mu0") = m_mu0, Rcpp::Named("m_mu") = m_mu, Rcpp::Named("Sigma_Mu1") = Sigma_Mu1, Rcpp::Named("Sigma_Mu0") = Sigma_Mu0, Rcpp::Named("Sigma_Mu") = Sigma_Mu);
  
}

// [[Rcpp::export]]
arma::vec FitTauStar(const arma::vec& m_mu, const arma::mat& Sigma_Mu, const arma::vec& nu_est, const arma::vec& VarNuStar, const double Atau_star, const double Btau_star, const double delta){
  
  const int p = m_mu.n_elem;
  
  arma::vec retVec(2);
  retVec(0) = Atau_star + 0.5*p;
  
  arma::vec Emu_sq = pow(m_mu,2) + Sigma_Mu.col(0);
  arma::vec Emumu = m_mu.rows(0,p-2)%m_mu.rows(1,p-1) + Sigma_Mu.submat(0,1,p-2,1);
  
  arma::vec expnu = exp(nu_est + 0.5*VarNuStar);
  const arma::vec b1_sq = 0.5*expnu/delta + 0.5*delta*exp(-nu_est + 0.5*VarNuStar) - 1.0;
  const arma::vec b0_sq = 0.5*expnu/delta;
  const arma::vec b0b1 = 0.5*ones(p-1) - b0_sq;
  
  retVec(1)  = Btau_star + 0.5*Emu_sq(0);
  for(int j = 0; j < (p-1); j++){
    
    retVec(1) += 0.5*(b1_sq(j)*Emu_sq(j)+2*b0b1(j)*Emumu(j)+b0_sq(j)*Emu_sq(j+1));
    
  }
  
  return retVec;
  
  
}

// [[Rcpp::export]]
List FitZ(const arma::mat& X, const arma::vec& y, const arma::vec& m_mu1, const arma::vec& m_mu0, const arma::vec& m_mu, const arma::vec& vw, const arma::mat& mEpsilonInvMAT, const arma::vec& zeta, const arma::vec& mS, const arma::mat& Sigma_S, const double& m_tau_inv, const double& delta){

  const int p = X.n_cols;
  const int n = y.n_elem;
  
  //QZ: col 0 = main-diag; col 1 = lower off-diag
  arma::cube QZ(p,2,n);
  arma::cube SigmaZ(p,2,n);
  QZ.zeros();
  arma::mat mZ(p,n);
  const arma::vec onesp = ones(p);
  arma::vec b1_sq(p-1);
  arma::vec b0_sq(p-1);
  arma::vec b1b0(p-1);
  arma::vec diagQ(p);
  diagQ.zeros();
  arma::vec RHS(p);
  const arma::vec b0_sq_base = exp(mS + 0.5*Sigma_S.col(0))/(2*delta);
  const arma::vec b1_sq_base = 0.5*delta*exp(-mS + 0.5*Sigma_S.col(0));
  
  for(int i = 0; i < n; i++){

    b0_sq = exp(zeta(i))*b0_sq_base;
    b1_sq = b0_sq + 0.5*delta*b1_sq_base*exp(-zeta(i)) - onesp.rows(0,p-2);
    b1b0 = 0.5*onesp.rows(0,p-2) - b0_sq;
    diagQ.zeros();
    diagQ(0) = 1.0;
    diagQ.rows(0,p-2) =  diagQ.rows(0,p-2) + b1_sq;
    diagQ.rows(1,p-1) =  diagQ.rows(1,p-1) + b0_sq;
    
    if(y(i)==1){
      
      QZ.subcube(0,0,i,p-1,0,i) = vw%(mEpsilonInvMAT.col(0)) + (onesp - vw)%mEpsilonInvMAT.col(2) + m_tau_inv*diagQ;
      RHS = vw%mEpsilonInvMAT.col(0)%(X.row(i).t()-m_mu1)  +  ((onesp-vw)%mEpsilonInvMAT.col(2))%((X.row(i)).t() - m_mu);
      
    }
    else{
      
      QZ.subcube(0,0,i,p-1,0,i) = vw%(mEpsilonInvMAT.col(1)) + (onesp - vw)%mEpsilonInvMAT.col(2) + m_tau_inv*diagQ;
      RHS = vw%mEpsilonInvMAT.col(1)%(X.row(i).t()-m_mu0)  +  ((onesp-vw)%mEpsilonInvMAT.col(2))%((X.row(i)).t() - m_mu);
      
    }
    
    
    
    QZ.subcube(0,1,i,p-2,1,i) = m_tau_inv*b1b0;
    mZ.col(i) = ThomasAlgo2(QZ.subcube(0,1,i,p-2,1,i), QZ.subcube(0,0,i,p-1,0,i), QZ.subcube(0,1,i,p-2,1,i), RHS);
    SigmaZ.slice(i) = BandedCholToInvNonStat2(QZ.subcube(0,0,i,p-1,0,i), QZ.subcube(0,1,i,p-1,1,i));
    
  }
  
  mZ = mZ.t();
  
  return Rcpp::List::create(Rcpp::Named("mZ") = mZ, Rcpp::Named("SigmaZ") = SigmaZ);
  
}

// [[Rcpp::export]]
List FitTauTS2(const arma::mat& mZ, const arma::cube& SigmaZ, const arma::vec& mS, const arma::mat& Sigma_S, const arma::vec& zeta, const double& atauTS, const double& btauTS, const double& delta){

  const int n = mZ.n_rows;
  const int p = mZ.n_cols;
  
  double coef1 = 0.0;
  double coef2 = 0.0;
  double coef3 = 0.0;
  double coef4 = 0.0;
  
  arma::vec tempEZsq(p);
  arma::mat tempSigmaZ(p,p);
  arma::vec Eb1sq(p-1);
  arma::vec Eb0sq(p-1);
  arma::vec Eb1b0(p-1);
  const arma::vec b0_sq_base = exp(mS + 0.5*Sigma_S.col(0))/(2*delta);
  const arma::vec b1_sq_base = 0.5*delta*exp(-mS + 0.5*Sigma_S.col(0));
  const arma::vec onesp = ones(p-1);
  
  for(int i = 0; i < n; i++){
    
    
    tempSigmaZ = SigmaZ.slice(i);
    tempEZsq = pow((mZ.row(i)).t(),2.0) + tempSigmaZ.col(0);
    Eb0sq = exp(zeta(i))*b0_sq_base;
    Eb1sq = Eb0sq + 0.5*delta*b1_sq_base*exp(-zeta(i)) - onesp;
    Eb1b0 = 0.5*onesp - Eb0sq;
    
    coef1 += tempEZsq(0);
    coef2 += accu(tempEZsq.rows(0,p-2)%Eb1sq);
    coef3 += accu(((mZ.submat(i,0,i,p-2)%mZ.submat(i,1,i,p-1)).t() + tempSigmaZ.submat(0,1,p-2,1))%Eb1b0);
    coef4 += accu(tempEZsq.rows(1,p-1)%Eb0sq);
    
  }
  
  double rtauTS = atauTS + 0.5*n*p;
  double stauTS = btauTS + 0.5*( coef1 + coef2 + 2*coef3 + coef4 );
  
  
  return Rcpp::List::create(Rcpp::Named("rtauTS") = rtauTS, Rcpp::Named("stauTS") = stauTS);
  
  
}

// [[Rcpp::export]]
List FitGamma(const arma::mat& X, const arma::vec& y, const arma::mat& mZ, const arma::cube& SigmaZ,  const arma::vec& m_mu1, const arma::vec& m_mu0, const arma::vec& m_mu, const arma::mat& Sigma_Mu1, const arma::mat& Sigma_Mu0, const arma::mat& Sigma_Mu, const arma::vec& vwOLD, const arma::mat& mEpsilonInvMAT, const arma::mat& mLogEpsilonMAT, const double& delta, const double& alpha, const double& beta){

  const int n = X.n_rows;
  const int p = X.n_cols;
  const int n1 = accu(y);
  const int n0 = n - n1;
  
  const arma::vec u = 0.5*n1*mLogEpsilonMAT.col(0) + 0.5*n0*mLogEpsilonMAT.col(1) - 0.5*n*mLogEpsilonMAT.col(2);
  arma::vec gColTotal(p);
  gColTotal.zeros();
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      
      if(y(i)==1){
        
        gColTotal(j) +=  mEpsilonInvMAT(j,0)*(pow((X(i,j) - m_mu1(j) - mZ(i,j)), 2) + Sigma_Mu1(j,0) + SigmaZ(j,0,i) ) - mEpsilonInvMAT(j,2)*(pow((X(i,j) - m_mu(j) - mZ(i,j)), 2) + Sigma_Mu(j,0) + SigmaZ(j,0,i) );
        
      }
      else{
        
        gColTotal(j) +=  mEpsilonInvMAT(j,1)*(pow((X(i,j) - m_mu0(j) - mZ(i,j)), 2) + Sigma_Mu0(j,0) + SigmaZ(j,0,i) ) - mEpsilonInvMAT(j,2)*(pow((X(i,j) - m_mu(j) - mZ(i,j)), 2) + Sigma_Mu(j,0) + SigmaZ(j,0,i) );
        
      }

    }
  }
  
  arma::vec logitvw = -u - 0.5*gColTotal;
  logitvw.rows(1,p-2) = logitvw.rows(1,p-2) - alpha*ones(p-2) + beta*(vwOLD.rows(0,p-3) + vwOLD.rows(2,p-1)); 
  logitvw(0) += -alpha + beta*vwOLD(1);
  logitvw(p-1) += -alpha + beta*vwOLD(p-2);
  const arma::vec vw = expit(logitvw);

  return Rcpp::List::create(Rcpp::Named("vw") = vw, Rcpp::Named("gColTotal") = gColTotal, Rcpp::Named("u") = u, Rcpp::Named("logitvw") = logitvw);
}

// [[Rcpp::export]]
arma::vec GradFuncPrecSVBCpp(const arma::vec& Ezsq, const arma::vec& Ezizi, const arma::vec v, const double& m_tau_inv, const double& delta, const double& tau2, const double& lambda){
  
  const int p = Ezsq.n_elem;
  
  const arma::vec db0sqdv = exp(v)/(2*delta);
  const arma::vec db1sqdv = db0sqdv - 0.5*delta*exp(-v);
  const arma::vec db0b1dv = -db0sqdv;
  
  const arma::vec LikelihoodPart = 0.5*ones(p-1) - 0.5*m_tau_inv*( db1sqdv%Ezsq.rows(0,p-2) + 2*db0b1dv%Ezizi + db0sqdv%Ezsq.rows(1,p-1));
  
  const double c0sq = lambda/(2.0*delta);
  const double c1sq = c0sq+ 0.5*delta/lambda - 1.0;
  const double c0c1 = 0.5*(1.0 - lambda/delta);
  
  arma::vec PriorPart(p-1,fill::zeros);
  PriorPart(0) += (1/tau2)*((1+c1sq)*v(0) + c0c1*v(1));
  PriorPart.rows(1,p-3) = (1/tau2)*((c0sq + c1sq)*v.rows(1,p-3) + c0c1*( v.rows(0,p-4) + v.rows(2,p-2) ));
  PriorPart(p-2) += (1/tau2)*(c0c1*v(p-3) +c0sq*v(p-2));
  
  const arma::vec ans = LikelihoodPart - PriorPart;
  return ans;
  
}

// [[Rcpp::export]]
double ELBONu(const arma::vec& Ezsq, const arma::vec& Ezizi, const arma::vec v, const double& m_tau_inv, const double& delta, const double& tau2, const double& lambda, const double& detOmega, const double& normS){

  const int p = Ezsq.n_elem;
  
  const arma::vec Onesvec(p-1,fill::ones);
  const arma::vec b0sq = exp(v)/(2*delta);
  const arma::vec b1sq = b0sq + 0.5*delta*exp(-v) - Onesvec;
  const arma::vec b0b1 = 0.5*Onesvec - b0sq;
  
  const double c0sq = lambda/(2.0*delta);
  const double c1sq = c0sq+ 0.5*delta/lambda - 1.0;
  const double c0c1 = 0.5*(1.0 - lambda/delta);
  
  const double LikelihoodPart = 0.5*accu(v) - 0.5*m_tau_inv*accu( b1sq%Ezsq.rows(0,p-2) + 2*b0b1%Ezizi + b0sq%Ezsq.rows(1,p-1) );
  arma::vec PriorPart(p-1,fill::zeros);
  PriorPart(0) += (1/tau2)*((1+c1sq)*v(0) + c0c1*v(1));
  PriorPart.rows(1,p-3) = (1/tau2)*((c0sq + c1sq)*v.rows(1,p-3) + c0c1*( v.rows(0,p-4) + v.rows(2,p-2) ));
  PriorPart(p-2) += (1/tau2)*(c0c1*v(p-3) +c0sq*v(p-2));
  const double PriorNumber = 0.5*accu(v%PriorPart);
  
  const double ans = LikelihoodPart - PriorNumber + 0.5*(p-1)*log(2.0*3.1416) - detOmega + 0.5*normS;
  return ans;
  
}

// [[Rcpp::export]]
List PrecSVBforNuCpp(arma::vec mZi, arma::mat SigmaZi, const double& m_tau_invTS, const double& delta, const double& tau2, const double& lambda, const arma::vec& mStart, const int B){

  const int p = SigmaZi.n_rows;
  const double rho = 0.7;
  const double epsilon = pow(10, -6);
  
  arma::vec mStore(p-1);
  arma::mat OmegaStore(p-1,2);
  arma::mat OmegaPrimeStore(p-1,2);
  
  arma::vec EgmSq(p-1, fill::zeros);
  arma::vec EDeltamSq(p-1, fill::zeros);
  
  arma::mat EgOmegaPrimeSq(p-1, 2, fill::zeros);
  arma::mat EDeltaOmegaPrimeSq(p-1, 2, fill::zeros);
  
  arma::vec ELBOstore(B);
  
  arma::vec mCurr = mStart;
  arma::mat OmegaCurr(p-1,2,fill::zeros);
  OmegaCurr.col(0) = ones(p-1);
  arma::mat OmegaPrimeCurr(p-1,2,fill::zeros);
  arma::vec vCurr(p-1,fill::zeros);
  arma::vec OmegasCurr(p-1,fill::zeros);
  
  arma::vec EziSq = pow(mZi, 2.0) + SigmaZi.col(0);
  arma::vec Ezizi = mZi.rows(0,p-2)%mZi.rows(1,p-1) + SigmaZi.submat(0,1,p-2,1);
  
  arma::mat gOmegaprime(p-1,2,fill::zeros);
  arma::mat DeltaOmegaPrimeCurr(p-1,2,fill::zeros);
  arma_rng::set_seed_random();
  arma::mat normsamp(p-1,B,fill::randn);
  
  arma::vec xi(p-1);
  const arma::vec zeroVec(p-2, fill::zeros);
  arma::vec epsVec(p-1);
  epsVec.fill(epsilon);
  arma::mat epsMat(p-1,2);
  epsMat.fill(epsilon);
  epsMat(p-2,1) = 0.0;
  arma::vec gradCurr(p-1);
  arma::vec gmCurr(p-1);
  arma::vec DeltamCurr(p-1);
  arma::vec OmegainvgmuCurr(p-1);
  arma::vec OmegainvsCurr(p-1);
  
  for(int iter = 0; iter < B; iter++){
    
    xi = ThomasAlgo2(zeroVec, OmegaCurr.col(0), OmegaCurr.submat(0,1,p-3,1), normsamp.col(iter));
    vCurr = mCurr + xi;
    gradCurr = GradFuncPrecSVBCpp(EziSq, Ezizi, vCurr, m_tau_invTS, delta, tau2, lambda);
    
    OmegasCurr.rows(1,p-2) = OmegaCurr.submat(1,0,p-2,0)%normsamp.submat(1,iter,p-2,iter) + OmegaCurr.submat(0,1,p-3,1)%normsamp.submat(0,iter,p-3,iter);
    OmegasCurr(0) = OmegaCurr(0,0)*normsamp(0,iter);
    gmCurr = gradCurr + OmegasCurr;
    EgmSq = rho*EgmSq + (1-rho)*pow(gmCurr,2.0);
    DeltamCurr = (sqrt(EDeltamSq + epsVec)/sqrt(EgmSq + epsVec))%gmCurr;
    EDeltamSq = rho*EDeltamSq + (1-rho)*pow(DeltamCurr,2.0);
    mCurr = mCurr + DeltamCurr;
    
    OmegainvgmuCurr = ThomasAlgo2(OmegaCurr.submat(0,1,p-3,1), OmegaCurr.col(0), zeroVec, gmCurr);
    gOmegaprime.col(0) = - xi%OmegainvgmuCurr;
    gOmegaprime.submat(0,1,p-3,1) = - xi.rows(1,p-2)%OmegainvgmuCurr.rows(0,p-3);
    gOmegaprime(p-2,1) = 0.0;
    gOmegaprime.col(0) = gOmegaprime.col(0)%OmegaCurr.col(0);
    EgOmegaPrimeSq = rho*EgOmegaPrimeSq + (1-rho)*pow(gOmegaprime,2.0);
    DeltaOmegaPrimeCurr = ((sqrt(EDeltaOmegaPrimeSq + epsMat))/(sqrt(EgOmegaPrimeSq + epsMat)))%gOmegaprime;
    DeltaOmegaPrimeCurr(p-2,1) = 0.0;
    EDeltaOmegaPrimeSq = rho*EDeltaOmegaPrimeSq + (1-rho)*pow(DeltaOmegaPrimeCurr, 2.0);
    OmegaPrimeCurr = OmegaPrimeCurr + DeltaOmegaPrimeCurr;
    OmegaCurr = OmegaPrimeCurr;
    OmegaCurr.col(0) = exp(OmegaPrimeCurr.col(0));
    
    
    ELBOstore(iter) = ELBONu(EziSq, Ezizi, vCurr, m_tau_invTS, delta, tau2, lambda, accu(OmegaPrimeCurr.col(0)), accu(pow(normsamp.col(iter),2.0)));
    
  }
  
  mStore = mCurr;
  OmegaStore = OmegaCurr;
  OmegaPrimeStore = OmegaPrimeCurr;
  
  return Rcpp::List::create(Rcpp::Named("mStore") = mStore, Rcpp::Named("OmegaStore") = OmegaStore, Rcpp::Named("OmegaPrimeStore") = OmegaPrimeStore, Rcpp::Named("ELBOStore") = ELBOstore);
  
}

// [[Rcpp::export]]
arma::vec GradFuncPrecScpp(const arma::vec& G1,  const arma::vec& G2, const arma::vec& G3, const double& c0sq, const double& c1sq, const double& c0c1, const double& m_tau_inv, const arma::vec& muS, const double& tau2, const arma::vec& Scurr, const double& delta, const arma::vec& zeta){

  const int p = G1.n_elem;
  const int n = zeta.n_elem;
  
  const arma::vec Stilde = Scurr - muS;
  arma::vec LikelihoodPart = 0.5*n*ones(p-1) - 0.5*m_tau_inv*( exp(Scurr)%G1.rows(0,p-2)/(2*delta) - 0.5*delta*exp(-Scurr)%G2.rows(0,p-2) - 2*exp(Scurr)%(G3/(2*delta)) + exp(Scurr)%G1.rows(1,p-1)/(2*delta)  );
  arma::vec PriorPart(p-1);
  PriorPart(0) = (1/tau2)*((1+c1sq)*Stilde(0) + c0c1*Stilde(1));
  PriorPart.rows(1,p-3) = (1/tau2)*((c1sq+c0sq)*Stilde.rows(1,p-3) + c0c1*(Stilde.rows(0,p-4) + Stilde.rows(2,p-2)) );
  PriorPart(p-2) = (1/tau2)*(c0c1*Stilde(p-3) + c0sq*Stilde(p-2));
  
  const arma::vec gradAns = LikelihoodPart - PriorPart;
  
  return gradAns;
  
}

// [[Rcpp::export]]
double ELBOSCpp(const arma::vec& G1,  const arma::vec& G2, const arma::vec& G3, const arma::vec& G4, const arma::vec& G5, const double& c0sq, const double& c1sq, const double& c0c1, const double& m_tau_inv, const arma::vec& muS, const double& tau2, const arma::vec& Scurr, const arma::vec& omegaCurrDiag, const arma::vec& normSamp, const double& delta, const arma::vec& zeta){
  
  const int p = G1.n_elem;
  const int nsamp = zeta.n_elem;
  const arma::vec Stilde = Scurr - muS;
  
  double logh = 0.5*nsamp*accu(Scurr) - 0.5*m_tau_inv*accu( exp(Scurr)%G1.rows(0,p-2)/(2*delta) + 0.5*delta*exp(-Scurr)%G2.rows(0,p-2) - G4.rows(0,p-2) + (G5 - exp(Scurr)%G3/(delta)) + exp(Scurr)%G1.rows(1,p-1)/(2*delta) );
  logh -= 0.5/tau2*( pow(Stilde(0),2) + c1sq*accu(pow(Stilde.rows(0,p-3),2)) + 2*c0c1*accu( Stilde.rows(0,p-3)%Stilde.rows(1,p-2) ) + c0sq*accu(pow(Stilde.rows(1,p-2),2)) );
  
  double ans = logh + 0.5*(p-1)*log(2*3.141593) - accu(log(omegaCurrDiag)) + 0.5*accu(pow(normSamp,2));
  
  return ans;
}

// [[Rcpp::export]]
List PrecSVBforSCpp(const arma::vec& G1,  const arma::vec& G2, const arma::vec& G3, const arma::vec& G4, const arma::vec& G5, const double& m_tau_inv, const arma::vec& muS, const double& tau2, const double& lambda, const arma::vec& mStart, const double& delta, const int& B, const arma::vec& zeta){

  const int p = G1.n_elem;
  
  const double rho =  0.7;
  const double epsilon = pow(10,-6);
  
  const double c0sq = lambda/(2*delta);
  const double c1sq = c0sq + delta/(2*lambda) - 1.0;
  const double c0c1 = 0.5 - c0sq;
  
  arma::vec mStore(p-1);
  arma::mat OmegaStore(p-1,2);
  arma::mat OmegaPrimeStore(p-1,2);
  
  arma::vec EgmSq(p-1, fill::zeros);
  arma::vec EDeltamSq(p-1, fill::zeros);
  
  arma::mat EgOmegaPrimeSq(p-1, 2, fill::zeros);
  arma::mat EDeltaOmegaPrimeSq(p-1, 2, fill::zeros);
  
  arma::vec ELBOstore(B);
  
  arma::vec mCurr = mStart;
  arma::mat OmegaCurr(p-1,2,fill::zeros);
  OmegaCurr.col(0) = ones(p-1);
  arma::mat OmegaPrimeCurr(p-1,2,fill::zeros);
  arma::vec SCurr(p-1,fill::zeros);
  arma::vec OmegasCurr(p-1,fill::zeros);
  
  arma::mat gOmegaprime(p-1,2,fill::zeros);
  arma::mat DeltaOmegaPrimeCurr(p-1,2,fill::zeros);
  arma_rng::set_seed_random();
  arma::mat normDraws(p-1,B,fill::randn);
  arma::vec xi(p-1);
  arma::vec normsamp(p-1);
  const arma::vec zeroVec(p-2, fill::zeros);
  arma::vec epsVec(p-1);
  epsVec.fill(epsilon);
  arma::mat epsMat(p-1,2);
  epsMat.fill(epsilon);
  epsMat(p-2,1) = 0.0;
  arma::vec gradCurr(p-1);
  arma::vec gmCurr(p-1);
  arma::vec DeltamCurr(p-1);
  arma::vec OmegainvgmuCurr(p-1);
  arma::vec OmegainvsCurr(p-1);
  
  for(int iter = 0; iter < B; iter++){
    
    normsamp = normDraws.col(iter);
    xi = ThomasAlgo2(zeroVec,OmegaCurr.col(0), OmegaCurr.submat(0,1,p-3,1), normsamp);
    SCurr = mCurr + xi;
    gradCurr = GradFuncPrecScpp(G1, G2, G3, c0sq, c1sq, c0c1, m_tau_inv, muS, tau2, SCurr, delta, zeta);
    OmegasCurr.rows(1,p-2) = OmegaCurr.submat(1,0,p-2,0)%normsamp.rows(1,p-2) + OmegaCurr.submat(0,1,p-3,1)%normsamp.rows(0,p-3);
    OmegasCurr(0) = OmegaCurr(0,0)*normsamp(0);
    gmCurr = gradCurr + OmegasCurr;
    EgmSq = rho*EgmSq + (1-rho)*pow(gmCurr,2.0);
    DeltamCurr = (sqrt(EDeltamSq + epsVec)/sqrt(EgmSq + epsVec))%gmCurr;
    EDeltamSq = rho*EDeltamSq + (1-rho)*pow(DeltamCurr,2.0);
    mCurr = mCurr + DeltamCurr;
    OmegainvgmuCurr = ThomasAlgo2(OmegaCurr.submat(0,1,p-3,1), OmegaCurr.col(0), zeroVec, gmCurr);
    gOmegaprime.col(0) = - xi%OmegainvgmuCurr;
    gOmegaprime.submat(0,1,p-3,1) = - xi.rows(1,p-2)%OmegainvgmuCurr.rows(0,p-3);
    gOmegaprime(p-2,1) = 0.0;
    gOmegaprime.col(0) = gOmegaprime.col(0)%OmegaCurr.col(0);
    EgOmegaPrimeSq = rho*EgOmegaPrimeSq + (1-rho)*pow(gOmegaprime,2.0);
    DeltaOmegaPrimeCurr = ((sqrt(EDeltaOmegaPrimeSq + epsMat))/(sqrt(EgOmegaPrimeSq + epsMat)))%gOmegaprime;
    DeltaOmegaPrimeCurr(p-2,1) = 0.0;
    EDeltaOmegaPrimeSq = rho*EDeltaOmegaPrimeSq + (1-rho)*pow(DeltaOmegaPrimeCurr, 2.0);
    OmegaPrimeCurr = OmegaPrimeCurr + DeltaOmegaPrimeCurr;
    OmegaCurr = OmegaPrimeCurr;
    OmegaCurr.col(0) = exp(OmegaPrimeCurr.col(0));
    ELBOstore(iter) = ELBOSCpp(G1, G2, G3, G4, G5, c0sq, c1sq, c0c1, m_tau_inv, muS, tau2, SCurr, OmegaCurr.col(0), normsamp, delta, zeta);
    
  }
  mStore = mCurr;
  OmegaStore = OmegaCurr;
  OmegaPrimeStore = OmegaPrimeCurr;
  
  return Rcpp::List::create(Rcpp::Named("mStore") = mStore, Rcpp::Named("OmegaStore") = OmegaStore, Rcpp::Named("OmegaPrimeStore") = OmegaPrimeStore, Rcpp::Named("ELBOStore") = ELBOstore);
  
}


// [[Rcpp::export]]
List FitPredZ(const arma::mat& X, const arma::vec& xi, const arma::vec& m_mu1, const arma::vec& m_mu0, const arma::vec& m_mu, const arma::vec& vw, const arma::mat& mEpsilonInvMAT, const arma::vec& zeta, const arma::vec& mS, const arma::mat& Sigma_S, const double& m_tau_inv, const double& delta){

  const int p = X.n_cols;
  const int n = xi.n_elem;
  
  //QZ: col 0 = main-diag; col 1 = lower off-diag
  arma::cube QZ(p,2,n);
  arma::cube SigmaZ(p,2,n);
  QZ.zeros();
  arma::mat mZ(p,n);
  const arma::vec onesp = ones(p);
  arma::vec b1_sq(p-1);
  arma::vec b0_sq(p-1);
  arma::vec b1b0(p-1);
  arma::vec diagQ(p);
  diagQ.zeros();
  arma::vec RHS(p);
  const arma::vec b0_sq_base = exp(mS + 0.5*Sigma_S.col(0))/(2*delta);
  const arma::vec b1_sq_base = 0.5*delta*exp(-mS + 0.5*Sigma_S.col(0));
  
  for(int i = 0; i < n; i++){
    
    b0_sq = exp(zeta(i))*b0_sq_base;
    b1_sq = b0_sq + 0.5*delta*b1_sq_base*exp(-zeta(i)) - onesp.rows(0,p-2);
    b1b0 = 0.5*onesp.rows(0,p-2) - b0_sq;
    diagQ.zeros();
    diagQ(0) = 1.0;
    diagQ.rows(0,p-2) =  diagQ.rows(0,p-2) + b1_sq;
    diagQ.rows(1,p-1) =  diagQ.rows(1,p-1) + b0_sq;
    
      QZ.subcube(0,0,i,p-1,0,i) = vw%(xi(i)*mEpsilonInvMAT.col(0) + (1-xi(i))*mEpsilonInvMAT.col(1)) + (onesp - vw)%mEpsilonInvMAT.col(2) + m_tau_inv*diagQ;
      RHS = (xi(i)*vw%mEpsilonInvMAT.col(0))%(X.row(i).t()-m_mu1) + ((1-xi(i))*vw%mEpsilonInvMAT.col(1))%(X.row(i).t()-m_mu0) +  ((onesp-vw)%mEpsilonInvMAT.col(2))%((X.row(i)).t() - m_mu);
    
    QZ.subcube(0,1,i,p-2,1,i) = m_tau_inv*b1b0;
    mZ.col(i) = ThomasAlgo2(QZ.subcube(0,1,i,p-2,1,i), QZ.subcube(0,0,i,p-1,0,i), QZ.subcube(0,1,i,p-2,1,i), RHS);
    SigmaZ.slice(i) = BandedCholToInvNonStat2(QZ.subcube(0,0,i,p-1,0,i), QZ.subcube(0,1,i,p-1,1,i));

  }
  
  mZ = mZ.t();
  
  return Rcpp::List::create(Rcpp::Named("mZ") = mZ, Rcpp::Named("SigmaZ") = SigmaZ);
  
}

// [[Rcpp::export]]
arma::vec FitPredY(const arma::mat& X, const arma::mat& mZ, const arma::cube& SigmaZ, const arma::vec& m_mu1, const arma::vec& m_mu0, const arma::mat& Sigma_Mu1, const arma::mat& Sigma_Mu0, const arma::mat& mEpsilonInvMAT, const arma::mat& mLogEpsilonMAT, const arma::vec& vw, const double& delta, const double& logRatio){

  const int m = X.n_rows;
  const int p = X.n_cols;
  
  const arma::vec wEpsilonInv1 = vw%mEpsilonInvMAT.col(0);
  const arma::vec wEpsilonInv0 = vw%mEpsilonInvMAT.col(1);
  const arma::vec DiffwEpsilonInv = wEpsilonInv1 - wEpsilonInv0;
  
  const double sumDiffLogEpsilon = accu(vw%(mLogEpsilonMAT.col(0) - mLogEpsilonMAT.col(1)));
  const double sumDiffEpsilonInvSigmaMu = accu( wEpsilonInv1%Sigma_Mu1.col(0) - wEpsilonInv0%Sigma_Mu0.col(0) );
  
  const arma::mat XminusZ = (X - mZ).t();
  double QDAtemp;
  double DiffSigmaZtemp;
  arma::vec logitxi(m);
  arma::vec VarZtemp(p);
  
  for(int i = 0; i < m; i++){
    
    VarZtemp = SigmaZ.subcube(0,0,i,p-1,0,i);
    QDAtemp = accu(wEpsilonInv1%pow( (XminusZ.col(i) - m_mu1), 2.0)) - accu(wEpsilonInv0%pow( (XminusZ.col(i) - m_mu0), 2.0));
    DiffSigmaZtemp = accu(DiffwEpsilonInv%VarZtemp);
    logitxi(i) = -0.5*sumDiffLogEpsilon - 0.5*QDAtemp - 0.5*sumDiffEpsilonInvSigmaMu - 0.5*DiffSigmaZtemp + logRatio;
    
  }
  
  const arma::vec xi = expit(logitxi);
  return xi;
}