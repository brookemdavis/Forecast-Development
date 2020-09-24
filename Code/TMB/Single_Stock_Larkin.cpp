#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_INTEGER(Priors);
  DATA_INTEGER(Scale);
  DATA_SCALAR(logA_mean);
  DATA_SCALAR(logA_sig);
  DATA_SCALAR(Sig_Gam_Dist);
  DATA_VECTOR(B_means);
  DATA_VECTOR(B_sigs);
  DATA_INTEGER(Bayes);
  
  PARAMETER(logA);
  PARAMETER_VECTOR(logB);
  PARAMETER(logSigma);
  
 // PARAMETER_VECTOR(logSgen);

  
  Type ans=0.0;
  Type sigma=exp(logSigma);
  int N_Obs = S.size();
  int N_Bs = logB.size();
  vector <Type> logR_Fit(N_Obs - N_Bs + 1);
  vector <Type> R_Fit(N_Obs - N_Bs + 1);
  vector <Type> logR_Pred(N_Obs - N_Bs + 1);
  Type A = exp(logA);
  vector <Type> B = exp(logB)/Scale;

  // vector <Type> SMSY( 4);  
   //vector <Type> logSMSY( 4);
  // vector <Type> Sgen = exp(logSgen);
  
  // for each stock need to start in year 4
  // indexing is not making sense
  for(int i=(N_Bs-1); i<N_Obs; i++){
      logR_Fit(i-N_Bs+1) = logA + log(S(i)) - exp(logB(0))*S(i) - 
                                        exp(logB(1))*S(i-1) - 
                                        exp(logB(2))*S(i-2) - 
                                        exp(logB(3))*S(i-3) ;
      ans += -dnorm(logR_Fit(i-N_Bs+1), logR(i),  sigma, true);
      R_Fit(i-N_Bs+1) = exp(logR_Fit(i-N_Bs+1))*Scale;
    }


  // Add priors
  if(Priors == 1){
    // Prior on logalpha -- no jacobian since same form as parameter
    ans -= dnorm(logA, logA_mean, logA_sig, true);
    // gamma prior on 1/sigma == inverse gamma on sigma, needs adjustment
    ans -= dgamma(pow(sigma, -2), Sig_Gam_Dist, 1/Sig_Gam_Dist, true);
    if(Bayes == 1){
      // Jacobian adjustment -- only needed if running as Bayesian model with tmbstan
      ans -= log(2)  - 2*logSigma;
    }
    // Follow jags model for now, normal dist on each beta
    for(int i=0; i<N_Bs; i++){
        ans -= dnorm(exp(logB(i)), B_means(i), B_sigs(i), true);
      if(Bayes == 1){
        // Jacobian adjustment -- only needed if running as Bayesian model with tmbstan
        ans -= logB(i);
      }
    }

  }


  SIMULATE {
    vector<Type> R_Pred(N_Obs-N_Bs+1);
    for(int i=0; i<(N_Obs-N_Bs+1); ++i){
      R_Pred(i) = exp(rnorm(logR_Fit(i), sigma)) * Scale;
    }
    REPORT(R_Pred);
  }
  
  //ADREPORT(SMSY);
  ADREPORT(A);
  //ADREPORT(Sgen);
  ADREPORT(sigma);
  ADREPORT(B);
  ADREPORT(R_Fit);
  REPORT(R_Fit);

  
  return ans;
  
}
