#include <TMB.hpp>
// Likelihood for a linear count model (model 1) 
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // Data time series and exogeneous covariate time series
  DATA_VECTOR(y);          
  DATA_VECTOR(x);        
  // Model parameters
  
  PARAMETER(d);
  PARAMETER(a); 
  PARAMETER(b); 
  PARAMETER(c);
  
  
  // Set initial values equal to stationary distribution
  Type ny0 = d/(1 - a - b);
  Type y0 = exp(ny0);
  // define 
  int n = y.size();
  vector<Type> ny(n);
  vector<Type> lambda(n);
  Type nll = 0.0;
  
  
  // ADREPORT on a, b, d
  ADREPORT(d);
  ADREPORT(a);
  ADREPORT(b);
  ADREPORT(c)
  
  // t=1
  ny(0) = d + a*ny0 + b*log(y0 + 1);
  lambda(0) = exp(ny(0));
  nll -= dpois(y(0), Type(lambda(0)), true);
  
  // t = 2,...,n
  for(int i = 1; i < n; i++){
    ny(i) = d + a*ny(i - 1) + b*log(y(i - 1) + 1) + c*x(i);
    lambda(i) = exp(ny(i));
    nll -= dpois(y(i), Type(lambda(i)), true);
  } //end of i
  return nll;
}