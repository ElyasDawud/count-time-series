#include <TMB.hpp>

// Likelihood for a linear count model (model 1) 
template<class Type>
  Type objective_function<Type>::operator() ()
{
  
  // Data
  DATA_VECTOR(y);          // timeseries vector
 
  
  // Parameters
  PARAMETER(d); 
  PARAMETER(a); 
  PARAMETER(b); 
   
  // ADREPORT on a, b, d
  ADREPORT(d);
  ADREPORT(a);
  ADREPORT(b);
  
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
  for(int i = 1; i < n; i++){
    lambda(i) = d + a*lambda(i - 1) + b*y(i - 1);
    // nll -= y(i - 1)*log(lambda(i - 1)) - lambda(i - 1);  
    nll -= dpois(y(i), lambda(i), true);
       } //end of i
  
  REPORT(lambda);
  
  return nll;
  }