//one_pops.stan

data {
int<lower=0> N; // number of data items
vector[N] eff; // eff vector
vector[N] rec; // rec vector
int<lower=0> N_new; //how many forecasts can be made?
vector[N_new] x_new;//whgat x values are there to make forecast can be made?
}
transformed data{
vector[N] x; // predictor vector
vector[N] y; // outcome vector

x=eff;
for(i in 1:N){
y[i]=log(rec[i]/eff[i]);
}

}
parameters {
real a_rick; // intercept
real b_rick; // coefficients for predictors
real<lower=0> sigma; // error scale
}
model {
y ~ normal(x * b_rick + a_rick, sigma); // likelihood
}
generated quantities {
vector[N_new] y_new;
vector[N_new] rec_pred;
for (n in 1:N_new){
y_new[n] = normal_rng(x_new[n] * b_rick + a_rick, sigma);
rec_pred[n]=x_new[n]*exp(y_new[n]);
}
}
