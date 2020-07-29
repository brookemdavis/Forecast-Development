//Single Stock Ricker Model

data {
int<lower=0> N; // number of data items
vector[N] R_Obs; // recruits
vector[N] S;   // spawners
real logA_mean;
real logA_sig;
real Sig_Gam_Dist;
real logSmax_mean;
real logSmax_sig;
}

parameters {
real logA; // intercept
real<lower=0> Smax; // coefficients for predictors
real<lower=0> sigma; // error scale
}



model {
vector[N] logR_Fit;
logR_Fit = logA + log(S) -  S/Smax;
R_Obs ~ lognormal(logR_Fit, sigma); // likelihood

//Priors
// lognormal prior on a
logA ~ normal( logA_mean, logA_sig);
// inverse gamma on variance (sigma^2), need jacobian adjustment
pow(sigma, 2) ~ inv_gamma(Sig_Gam_Dist, Sig_Gam_Dist);
// Jacobian adjustment
target += log(2*sigma); // log|x^2 d/dx|
// Lognormal prior on Smax (normal on log(Smax)) 
Smax ~ lognormal( logSmax_mean, logSmax_sig);
}

generated quantities {
vector[N] R_Pred;
vector[N] R_Fit;
for (i in 1:N){
 R_Pred[i] = lognormal_rng(A + B*log(S[i]), sigma);
 R_Fit[i] = exp(A + B*log(S[i]));
}
}
