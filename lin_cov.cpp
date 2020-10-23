#include <TMB.hpp>

// Likelihood for a linear count model (model 1) 
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // Data
  DATA_VECTOR(y);          // timeseries vector
  DATA_VECTOR(x); 
  
  
  // Parameters
  PARAMETER(d); 
  PARAMETER(a); 
  PARAMETER(b); 
  PARAMETER(c);
  
  // ADREPORT on a, b, d
  ADREPORT(d);
  ADREPORT(a);
  ADREPORT(b);
  ADREPORT(c)
  
  // Set initial values 
  Type y0 = d/(1 - a - b);
  Type lambda0 = y0;  
  int n = y.size();
  vector<Type> lambda(n);
  Type nll = 0.0;
  
  // t=1
  lambda(0) = d + a*lambda0 + b*y0;
  // nll -= y(0)*log(lambda(0)) - lambda(0);
  nll -= dpois(y(0), lambda(0), true);
  // t = 2,...,n
  for(int t = 1; t < n; t++){
    lambda(t) = d + a*lambda(t - 1) + b*y(t - 1) +c*x(t);
    // nll -= y(i - 1)*log(lambda(i - 1)) - lambda(i - 1);  
    nll -= dpois(y(t), lambda(t), true);
  } //end of i
  
  REPORT(lambda);
  
  return nll;
}