//Single Stock Ricker Model

data {
int<lower=0> N; // number of data items
vector[N] R_Obs; // recruits
vector[N] S;   // spawners
real A_mean;
real A_sig;
real Sig_Gam_Dist;
real B_mean;
real B_sig;
}

parameters {
real A; // intercept
real B; // coefficients for predictors
real<lower=0> sigma; // error scale
}



model {
vector[N] logR_Fit;
logR_Fit = A + B*log(S);
R_Obs ~ lognormal(logR_Fit, sigma); // likelihood

//Priors
// normal prior on a
A ~ normal( A_mean, A_sig);
// inverse gamma on variance (sigma^2), need jacobian adjustment
pow(sigma, 2) ~ inv_gamma(Sig_Gam_Dist, Sig_Gam_Dist);
// Jacobian adjustment
target += log(2*sigma); // log|x^2 d/dx|
// Normal prior on B 
B ~ normal( B_mean, B_sig);
}

generated quantities {
vector[N] R_Pred;
vector[N] R_Fit;
for (i in 1:N){
 R_Pred[i] = lognormal_rng(A + B*log(S[i]), sigma);
 R_Fit[i] = exp(A + B*log(S[i]));
 
}
}
