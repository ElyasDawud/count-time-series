#include <TMB.hpp>
// Likelihood for a bivariate poisson linear count model (model 1) 
template<class Type>
Type objective_function<Type>::operator() ()
{
  // bivariate posisson distributed count time series,Nx2 matrix, with a bivariate covariate t.s
  DATA_MATRIX(Y);
  DATA_MATRIX(X);
  // parameters
  PARAMETER_VECTOR(p);
  
  // split the data into two column vectors,
  vector<Type> Y1 = Y.col(0);    // 1st col of Y
  vector<Type> Y2 = Y.col(1);    // 2dn col of Y
  vector<Type> X1 = X.col(0);    // 1st col of X
  vector<Type> X2 = X.col(1);    // 2dn col of X
  
  
  int n = Y1.size();              // Length of the data
  
  
  vector<Type> loglik1(n);
  vector<Type> loglik2(n);
  vector<Type> loglik(n);
  
  // report parameter estimates with standard errors
  ADREPORT(p);
  
  vector<Type> nu1(n);
  vector<Type> nu2(n);
  vector<Type> lambda1(n);
  vector<Type> lambda2(n);
  //Some initial values
  loglik1(0) = 0;
  loglik2(0) = 0;
  loglik(0) = 0;
  nu1(0) = 0;
  nu2(0) = 0;
  lambda1(0) = 1;
  lambda2(0) = 1;
  for(int t=1;t<n;t++){
    nu1(t) = p(0) + p(2)*nu1(t-1) + p(4)*nu2(t-1) + p(6)*log(Y1(t-1)+1) + p(8)*log(Y2(t-1)+1)+p(10)*X1(t)+p(12)*X2(t);
    lambda1(t) = exp(nu1(t));
    nu2(t) = p(1) + p(3)*nu1(t-1) + p(5)*nu2(t-1) + p(7)*log(Y1(t-1)+1) + p(9)*log(Y2(t-1)+1)+p(11)*X1(t)+p(13)*X2(t);
    lambda2(t) = exp(nu2(t));
    if(lambda1(t)<=0){
      loglik1(t)=0;
    }
    else if(lambda1(t)>0){
      loglik1(t) = -Y1(t)*log(lambda1(t)) + lambda1(t);
    }
    if(lambda2(t)<=0){
      loglik2(t) = 0;
    }
    else if(lambda2(t)>0){
      loglik2(t) = -Y2(t)*log(lambda2(t)) + lambda2(t);
    }
    loglik(t) = loglik1(t) + loglik2(t);
  }
  return(loglik.sum());
}