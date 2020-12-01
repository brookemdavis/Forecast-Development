#include <TMB.hpp>

// Set up Lambert's W function to use to calculate SMSY
// Code taken from https://kaskr.github.io/adcomp/lambert_8cpp_source.html
// Step 1: Code up a plain C version
// Double version of Lambert W function
double LambertW(double x) {
  double logx = log(x);
  double y = (logx > 0 ? logx : 0);
  int niter = 100, i=0;
  for (; i < niter; i++) {
    if ( fabs( logx - log(y) - y) < 1e-9) break;
    y -= (y - exp(logx - y)) / (1 + y);
  }
  if (i == niter) Rf_warning("W: failed convergence");
  return y;
}

TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  LambertW
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = LambertW(tx[0]); // Call the 'double' version
,
// ATOMIC_REVERSE
Type W  = ty[0];                    // Function value from forward pass
Type DW = 1. / (exp(W) * (1. + W)); // Derivative
px[0] = DW * py[0];                 // Reverse mode chain rule
)
  
  // Scalar version
  template<class Type>
  Type LambertW(Type x){
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return LambertW(tx)[0];
  }
  
  
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logR);
  DATA_INTEGER(Priors);
  DATA_INTEGER(BiasCorr);
  DATA_INTEGER(Scale);
  DATA_SCALAR(logSmax_mean);
  DATA_SCALAR(logSmax_sig);
  DATA_SCALAR(logA_mean);
  DATA_SCALAR(logA_sig);
  DATA_SCALAR(Sig_Gam_Dist);
  DATA_INTEGER(Bayes);

  PARAMETER(logA);
  PARAMETER(logSigma);
  PARAMETER(logSmax);

  Type ans=0;
  Type sigma=exp(logSigma);
  int N_Obs = S.size(); 
  vector<Type> logR_Fit(N_Obs);
  vector<Type> R_Fit(N_Obs);
  Type Smax = exp(logSmax);
  Type B = 1/(Smax/Scale);
  Type A = exp(logA);
  
  
  for(int i=0; i<N_Obs; ++i){
    logR_Fit(i) = logA + log(S(i)) -  S(i)/Smax;
    if(BiasCorr == 0){
      ans += -dnorm(logR(i), logR_Fit(i), sigma, true);
    }
    if(BiasCorr == 1){
      // Back-transformation bias correction. These two versions give the same answer:
      //ans += -dnorm(logR(i)-logR_Fit(i), -pow(sigma,2)/Type(2), sigma, true);
      ans += -dnorm(logR(i) - logR_Fit(i) + pow(sigma,2)/Type(2), Type(0), sigma, true);
      
    }
    R_Fit(i) = exp(logR_Fit(i))*Scale;
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
   
  // Lognormal prior on Smax (normal on logSmax) 
   ans -= dnorm(logSmax, logSmax_mean, logSmax_sig, true);
 }
 // Calculate SMSY using Lambert's W function
 // Approach from Scheurell 2016
 Type SMSY =  (1 - LambertW(exp(1-logA)) ) / B ;
 

 SIMULATE {
   vector<Type> R_Pred(N_Obs);
   for(int i=0; i<N_Obs; ++i){
     R_Pred(i) = exp(rnorm(logR_Fit(i), sigma)) * Scale;
   }
   REPORT(R_Pred);
 }
  
  ADREPORT(A);
  ADREPORT(Smax);
  ADREPORT(sigma);
  ADREPORT(R_Fit);
  ADREPORT(SMSY);
  REPORT(R_Fit);
  return ans;
  
}
