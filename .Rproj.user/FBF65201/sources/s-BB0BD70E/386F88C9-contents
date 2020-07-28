#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SCALAR(Sig_Gam_Dist);
  DATA_SCALAR(A_mean);
  DATA_SCALAR(A_sig);

  // We want a N(0,1) prior on the standard deviation. Do it with
  // and without using external bounds. WIthout requires
  // exponentiating here and then adding a Jacobian
  // adjustment. That is done automatically when putting bounds
  // in R.  
  PARAMETER(sigma1);    // use external bounds
  PARAMETER(logsigma2); // no external bounds
  // Can also do Inverse Gamma
  PARAMETER(tau3);	// use external bounds
  PARAMETER(logsigma4);	// no external bounds
  PARAMETER(logA);

  Type sigma2=exp(logsigma2);
  Type sigma4=exp(logsigma4);
  Type var=pow(sigma4, Type(2));
  Type tau4=1/var;
  Type A = exp(logA);
  
  Type ans = 0;

  ans -= dnorm(sigma1, Type(0.0), Type(1.0), true);
  ans -= dnorm(sigma2, Type(0.0), Type(1.0), true);
  // Jacobian adjustment for sigma2
  ans -= logsigma2;
  
  // gamma distribution on tau == inverse gamma on variance
  // gamma in TMB is shape and scale (not rate, like in JAGS)
  ans -= dgamma(tau3, Sig_Gam_Dist, 1/Sig_Gam_Dist, true);
  ans -= dgamma(tau4, Sig_Gam_Dist, 1/Sig_Gam_Dist, true);
  // Jacobian adjustment
  ans -= log(2) - 2*logsigma4;
  
  
  // normal dist on A
  ans -= dnorm(A, A_mean, A_sig, true);
  ans -= logA;
 
  return ans;
}
