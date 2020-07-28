# Test Ricker, Larkin, Recursive, Sibling TMB models

library(TMB)
#library(tmbstan)
library(dplyr)
library(ggplot2)
source("Code/Functions.R")

#================================================================================

# Ricker model

# Simulate some basic Ricker data

SimData <- fake_SR_data(leng=50, age=4, Sig_Ricker = 0.2, true_a = rnorm(1, 5, 2), true_b=1/100000,
                          hr_min = 0.2, hr_max = 0.8)

SimDataDF <- data.frame(S = round(SimData$esc), R = (SimData$rec), Year = 1:length(SimData$esc))

ggplot(SimDataDF, aes(x=S, y=R)) + geom_point() + coord_fixed()

# only need to compile if changed model
#dyn.unload(dynlib("Code/TMB/Single_Stock_Ricker"))
TMB::compile("Code/TMB/Single_Stock_Ricker.cpp") # can't seem to compile if have TMBstan loaded
dyn.load(dynlib("Code/TMB/Single_Stock_Ricker"))

# want to scale down all data so close to 0
Scale <- 10^(floor(log(mean(SimDataDF$R), 10)))

data <- list()
data$S <- SimDataDF$S/Scale 
data$logR <- log(SimDataDF$R/Scale)
data$Priors <- 0
# need to have these in Data even though aren't used
data$logA_mean <- 0
data$logA_sig <- 0
data$Sig_Gam_Dist <- 0
data$Smax_mean <- 0
data$Smax_sig <- 0
data$Scale <- Scale


param <- list()
param$logA <- 1
param$logB <- log(1/(quantile(SimDataDF$S, 0.8)/Scale) )
param$logSigma <- -2

# Now Fit TMB model
obj <- MakeADFun(data, param, DLL="Single_Stock_Ricker", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))

# pull out estimates from ML fit
# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

# pull out predictions
R_Ests <- All_Ests[grepl("R_Pred", All_Ests$Param),  ] [, -3 ]
names(R_Ests) <- c("R_Pred", "StdErr")

library(tmbstan)
# Now run as MCMC? -- would only do this for final year
fitmcmc <- tmbstan(obj, chains=3, iter=100000, init=list(opt$par))
                    #control = list(adapt_delta = 0.89))

xx<-as.matrix(fitmcmc)


mc <- extract(fitmcmc, pars=names(obj$par),  
              inc_warmup=FALSE, permuted=FALSE)  

xx<-sampling(fitmcmc)

# extract fitted values -- this doesn't work, how to pull out predictions? -- maybe get posterior predictive?
Preds <- extract(fitmcmc, pars=All_Ests$Param[grepl("R_Pred", All_Ests$Param)],  
              inc_warmup=FALSE, permuted=FALSE) 

# not sure if I want to dwell on this too much, since we will only need posterior around one or two values,
# which can just calculate manually for each run?


# add to simDataDF
SimDataDF <- cbind(SimDataDF, R_Ests)
SimDataDF$CI_low <- SimDataDF$R_Pred - 1.96*SimDataDF$StdErr
SimDataDF$CI_up <- SimDataDF$R_Pred + 1.96*SimDataDF$StdErr
SimDataDF$True <- SimData$true_a * SimDataDF$S * exp( -SimData$true_b * SimDataDF$S )

ggplot(SimDataDF, aes(S, R, ymin = CI_low, ymax = CI_up)) +
  geom_point() +
  geom_line(aes(y=R_Pred)) +
  geom_ribbon( alpha = 0.1) +
  geom_line(aes(y=True), col = "red")
  theme_bw()
  
  
# want to add uncertainty bars -- then compare to version with priors, and JAGS version

# make it so model makes predictions so can get CI's around them -- can always remove later for forecast
   # will replace with forecast using y-4 and y-5



