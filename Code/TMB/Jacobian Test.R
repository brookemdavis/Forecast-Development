# write toy TMB model to check that have Jacobian right

library(TMB)
library(pscl) # use this package for nice inverse gamma functions
library(truncnorm)
dyn.unload(dynlib("Code/TMB/SigmaTest"))
TMB::compile("Code/TMB/SigmaTest.cpp")
dyn.load(dynlib("Code/TMB/SigmaTest"))

Sig_Gam_Dist <- 1
A_mean <- 5
A_sig <- 1

data <- list(Sig_Gam_Dist=Sig_Gam_Dist, A_mean = A_mean, A_sig = A_sig)
param <- list(sigma1=1, logsigma2=0, tau3=1, logsigma4=1, logA = 1)

obj <- MakeADFun(data, param, DLL="SigmaTest" )
obj$fn()
summary(sdreport(obj))

library(tmbstan)
library(shinystan)
fitmcmc <- tmbstan(obj, chains=3, iter=200000, thin=10,
                   lower=c(0,-Inf,0,-Inf, -Inf), upper=c(Inf, Inf, Inf,Inf, Inf))
# pull out posterior vals
post <- as.data.frame(fitmcmc)

par(mfrow=c(2,2))
x <- seq(1e-5,10, len=500)
bb <- 50
hist(post$sigma1, freq=FALSE, breaks=bb)
lines(x, dtruncnorm(x, 0,Inf))
hist(exp(post$logsigma2), freq=FALSE, breaks=bb)
lines(x, dtruncnorm(x, 0,Inf))

# plot his of tau, which should be gamma dist
hist(post$tau3, freq=FALSE, breaks=bb)
lines(x, dgamma(x,1,1))
hist(1/exp(post$logsigma4)^2, freq=FALSE, breaks=bb)
lines(x, dgamma(x, 1,1))

# plot hist of variance, which should be inverse gamma dist.
VarPost <- exp(2*post$logsigma4)
# variance posterior has some massively large values, in order to look
# at histogram need to cut out some of these large values
# 90% of the posterior density is below 10 -- use this as cutoff
hist(VarPost[VarPost<10], freq=FALSE, breaks=bb)
lines(x, densigamma(x, 1,1))

x <- seq(1, 10, len=500)
bb <- 50
hist(exp(post$logA), freq=FALSE, breaks=bb)
lines(x, dnorm(x, A_mean, A_sig))
# looks good

# simulate gamma on tau to look at dist and median
tau_sim <- rgamma(10000, shape = Sig_Gam_Dist, rate = Sig_Gam_Dist)
median(tau_sim)
var <- 1/tau_sim
median(var)
median(VarPost)
# looks good!

# simulate inverse gamma directly on variance
var_sim <- rigamma(100000, alpha = Sig_Gam_Dist, beta = Sig_Gam_Dist)
median(var_sim) # same as above
median(VarPost)
# looks good!






