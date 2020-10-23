#include <TMB.hpp>
// Likelihood for a bivariate poisson linear count model (model 1) 
template<class Type>
Type objective_function<Type>::operator() ()
{
  // bivariate posisson distributed count time series,Nx2 matrix
  DATA_MATRIX(Y);
  // parameters
  PARAMETER(d1);
  PARAMETER(d2);
  PARAMETER(a11);
  PARAMETER(a12);
  PARAMETER(a21);
  PARAMETER(a22);
  PARAMETER(b11);
  PARAMETER(b12);
  PARAMETER(b21);
  PARAMETER(b22);
  
  // split the data into two vectors
  vector<Type> Y1 = Y.col(0);    // 1st col of Y
  vector<Type> Y2 = Y.col(1);    // 2dn col of Y
  
  int n = Y1.size();              // Length of the data
  
  
  vector<Type> loglik1(n);
  vector<Type> loglik2(n);
  vector<Type> loglik(n);

  // report parameter estimates with standard errors
  ADREPORT(d1);
  ADREPORT(d2);
  ADREPORT(a11);
  ADREPORT(a12);
  ADREPORT(a21);
  ADREPORT(a22);
  ADREPORT(b11);
  ADREPORT(b12);
  ADREPORT(b21);
  ADREPORT(b22);
  
  vector<Type> lambda1(n);
  vector<Type> lambda2(n);
  lambda1(0) = 1;
  lambda2(0) = 1;
  loglik1(0) = 0;
  loglik2(0) = 0;
  loglik(0) = 0;
  
  for(int t=1;t<n;t++){
    lambda1(t) = d1 + a11*lambda1(t-1) + a12*lambda2(t-1) + b11*Y1(t-1) + b12*Y2(t-1);
    lambda2(t) = d2 + a21*lambda1(t-1) + a22*lambda2(t-1) + b21*Y1(t-1) + b22*Y2(t-1);
    if(lambda1(t)<=0){
      loglik1(t)=0;
    }
    else if(lambda1(t)>0){
      loglik1(t) = -Y1(t)*log(lambda1(t))+lambda1(t);
    }
    if(lambda2(t)<=0){
      loglik2(t) = 0;
    }
    else if(lambda2(t)>0){
      loglik2(t) = -Y2(t)*log(lambda2(t))+lambda2(t);
    }
    loglik(t) = loglik1(t) + loglik2(t);
  }
  return(loglik.sum());
}