#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SCALAR(Sig_Gam_Dist);
  
  PARAMETER(logSigma);

  Type ans=0;
  Type sigma=exp(logSigma);
  
 
  ans -= dgamma(pow(sigma, -2), Sig_Gam_Dist, 1/Sig_Gam_Dist, true);
  // Jacobian adjustment
  ans -= log(2)  - 3*logSigma;
 
  ADREPORT(sigma);
  
  return ans;
  
}
