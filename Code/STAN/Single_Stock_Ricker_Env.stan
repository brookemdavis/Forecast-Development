//Single Stock Ricker Model

data {
int<lower=0> N; // number of data items
vector[N] R_Obs; // recruits
vector[N] S;   // spawners
vector[N] env; // Environmental covariate
real logA_mean;
real logA_sig;
real Sig_Gam_Dist;
real logSmax_mean;
real logSmax_sig;
real g_mean;
real g_sig;
}

parameters {
real logA; // intercept
real<lower=0> Smax; // coefficients for predictors
real<lower=0> sigma; // error scale
real g; //env. covariate param.
}



model {
vector[N] logR_Fit;
logR_Fit = logA + log(S) -  S/Smax + g * env;
R_Obs ~ lognormal(logR_Fit, sigma); // likelihood

//Priors
// lognormal prior on a
logA ~ normal( logA_mean, logA_sig);
// normal prior on g
g ~ normal( g_mean, g_sig);
// inverse gamma on variance (sigma^2), need jacobian adjustment
pow(sigma, 2) ~ inv_gamma(Sig_Gam_Dist, Sig_Gam_Dist);
// Jacobian adjustment
target += log(2*sigma); // log|x^2 d/dx|
// Lognormal prior on Smax (normal on log(Smax)) 
Smax ~ lognormal( logSmax_mean, logSmax_sig);
}

generated quantities {
  real SMSY = logA*(0.5-0.07*logA)*Smax;
vector[N] R_Pred;
vector[N] R_Fit;
for (i in 1:N){
 R_Pred[i] = lognormal_rng(logA + log(S[i]) -  S[i]/Smax + g * env, sigma);
 R_Fit[i] = exp(logA + log(S[i]) -  S[i]/Smax + g * env);
}
}
