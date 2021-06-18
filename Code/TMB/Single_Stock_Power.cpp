#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_INTEGER(Priors);
  DATA_INTEGER(Scale);
  DATA_SCALAR(A_mean);
  DATA_SCALAR(A_sig);
  DATA_SCALAR(Sig_Gam_Dist);
  DATA_SCALAR(B_mean);
  DATA_SCALAR(B_sig);
  DATA_INTEGER(Bayes);
  
  PARAMETER(logA);
  PARAMETER(logSigma);
  PARAMETER(logB);
  
  Type ans=0;
  Type sigma=exp(logSigma);
  int N_Obs = S.size(); 
  vector<Type> logR_Fit(N_Obs);
  vector<Type> R_Fit(N_Obs);
  Type B = exp(logB);
  Type A = exp(logA);
  
  
  for(int i=0; i<N_Obs; ++i){
    logR_Fit(i) = A + B*log(S(i));
    ans += -dnorm(logR(i), logR_Fit(i), sigma, true);
    R_Fit(i) = exp(logR_Fit(i))*Scale;
  }
  
  // Add priors
  if(Priors == 1){
    // Normal prior on A -- this is form in current forecast, will leave as is for now
    ans -= dnorm(exp(logA), A_mean, A_sig, true);
    if(Bayes == 1){
      // need Jacobian since logA is parameter log | exp(x) d/dx | = x
      ans -= logA;
    }
    // gamma prior on 1/sigma == inverse gamma on sigma, needs adjustment 
    ans -= dgamma(pow(sigma, -2), Sig_Gam_Dist, 1/Sig_Gam_Dist, true);
    if(Bayes == 1){
      // Jacobian adjustment
      ans -= log(2)  - 2*logSigma;
    }
    // Normal prior on B
    ans -= dnorm(exp(logB), B_mean, B_sig, true);
    if(Bayes == 1){
      // need Jacobian since logB is parameter log | exp(x) d/dx | = x
      ans -= logB;
    }
  }
  
  SIMULATE {
    vector<Type> R_Pred(N_Obs);
    for(int i=0; i<N_Obs; ++i){
      R_Pred(i) = exp(rnorm(logR_Fit(i), sigma)) * Scale;
    }
    REPORT(R_Pred);
  }
  
  ADREPORT(A);
  ADREPORT(B);
  ADREPORT(sigma);
  ADREPORT(R_Fit);
  REPORT(R_Fit);
  return ans;
  
}